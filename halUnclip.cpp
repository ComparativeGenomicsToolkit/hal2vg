/*
 * Copyright (C) 2016 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

// Convert clipped sequences (like chr1_sub_110000_22220000) back to their original states

//#define debug

#include <cstdlib>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <deque>
#include <unordered_map>
#include <unordered_set>

#include "hal.h"
#include "commonC.h"
#include "bioioC.h"
#include "subpaths.h"

using namespace std;
using namespace hal;

static void initParser(CLParser* optionsParser) {
    optionsParser->addArgument("inFile", "input HAL file");
    optionsParser->addArgument("seqFile", "cactus-style seqfile. 1st col=genome name, 2nd col=(original) fasta file.  only local paths supported");
    optionsParser->addArgument("outFile", "output HAL file");
    optionsParser->addOptionFlag("progress",
                                 "show progress",
                                 false);
    optionsParser->addOptionFlag("validate",
                                 "run a (non-exhaustive) check on the output",
                                 false);
    optionsParser->addOption("targetGenomes",
                             "comma-separated (no spaces) list of target genomes "
                             "(others are not unclipped) (all leaves if empty)",
                             "\"\"");
    optionsParser->setDescription("Fill back clipped sequence (removed by cactus-preprocess) using the original fasta files"
        ". Star trees only");
}

static vector<string> split_delims(const string &s, const string& delims) {
    vector<string> elems;
    size_t start = string::npos;
    for (size_t i = 0; i < s.size(); ++i) {
        if (delims.find(s[i]) != string::npos) {
            if (start != string::npos && i > start) {
                elems.push_back(s.substr(start, i - start));
            }
            start = string::npos;
        } else if (start == string::npos) {
            start = i;
        }
    }
    if (start != string::npos && start < s.size()) {
        elems.push_back(s.substr(start, s.size() - start));
    }
    return elems;
}

// c++ wrapper for sonlib -- load fasta file into memory
static unordered_map<string, pair<size_t, string>> read_fasta(const string& fa_path) {
    FILE* fa_file = fopen(fa_path.c_str(), "r");
    if (!fa_file) {
        cerr << "Unable to open fastat file: " << fa_path << endl;
        exit(1);
    }

    List* seqs = constructEmptyList(0, free);
    List* seq_lens = constructEmptyList(0, free);
    List* seq_names = constructEmptyList(0, free);

    fastaRead(fa_file, seqs, seq_lens, seq_names);

    // todo: should be done once, but sonlib fasta reading so slow it odesn't matter
    vector<unsigned char> cmap(numeric_limits<unsigned char>::max());
    for (unsigned char i = 0; i < cmap.size(); ++i) {
        switch (i) {
        case 'a':
        case 'c':
        case 'g':
        case 't':
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            cmap[i] = i;
            break;
        default:
            cmap[i] = 'N';
            break;
        }
    }

    unordered_map<string, pair<size_t, string>> fa_info; 
    for (int64_t i = 0; i < seqs->length; ++i) {
        string name = (char*)seq_names->list[i];
        size_t len = (size_t)listGetInt(seq_lens, i);
        string seq = (char*)seqs->list[i];
        for (size_t j = 0; j < seq.length(); ++j) {
            // hal doesn't like non-acgtn characters
            seq[j] = cmap[seq[j]];
        }
        fa_info[name] = make_pair(len, seq);
    }

    destructList(seqs);
    destructList(seq_lens);
    destructList(seq_names);

    fclose(fa_file);

    return fa_info;
}

// do a pass over the seqfile to get the total lengths of every sequence
static unordered_map<string, size_t> get_dimensions_from_seqfile(const string& seqfile_path, const unordered_set<string>& target_set) {
    unordered_map<string, size_t> seq_map;
    
    ifstream seqfile(seqfile_path);
    if (!seqfile) {
        cerr << "[halUnclip]: Unable to open seqfile: " << seqfile_path << endl;
        exit(1);
    }

    string buffer;
    while (getline(seqfile, buffer)) {
        vector<string> toks = split_delims(buffer, " \t");
        if (toks.size() == 2) {
            string name = toks[0];
            string fa_path = toks[1];
            if (target_set.count(name)) {
                unordered_map<string, pair<size_t, string>> fa_info = read_fasta(fa_path);
                for (auto& fi : fa_info) {
                    seq_map[name + "." + fi.first] = fi.second.first;
                }
            }
        }
    }
    
    return seq_map;
}

static unordered_map<string, vector<Sequence::Info>> get_filled_dimensions(AlignmentConstPtr alignment, unordered_map<string, size_t>& seq_d,
                                                                           const unordered_set<string>& target_set, bool progress) {

    unordered_map<string, vector<Sequence::Info>> dim_map;

    // copy root exactly as is
    vector<Sequence::Info>& root_dims = dim_map[alignment->getRootName()];
    const Genome* root_genome = alignment->openGenome(alignment->getRootName());
    for (SequenceIteratorPtr seqIt = root_genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
        const Sequence *sequence = seqIt->getSequence();
        root_dims.push_back(Sequence::Info(sequence->getName(), sequence->getSequenceLength(), sequence->getNumTopSegments(), sequence->getNumBottomSegments()));
    }
    
    vector<string> names = alignment->getChildNames(alignment->getRootName());

    for (const string& name : names) {
        const Genome* genome = alignment->openGenome(name);
        vector<Sequence::Info>&  dimensions = dim_map[name];
        if (progress) {
            cerr << "[halUnclip]: Scanning dimensions of genome " << genome->getName() << endl;
        }

        // map base name to sequence fragments
        unordered_map<string, vector<const Sequence*>> frag_map;

        // pass 1, map all hal sequences back to their base name and check that they correspond to a fasta sequence
        for (SequenceIteratorPtr seqIt = genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
            const Sequence *sequence = seqIt->getSequence();
            string sequence_name = sequence->getName();
            string parsed_name = target_set.count(genome->getName()) ? parse_subpath_name(sequence_name) : sequence_name;
            string full_name = genome->getName() + "." + parsed_name;
            size_t fa_len = sequence->getSequenceLength();
            if (!seq_d.count(full_name)) {
                if (parsed_name != sequence_name) {
                    cerr << "[halUnclip]: Unable to find sequence (from HAL) " << full_name << " in dimension map from input fasta" << endl;
                    exit(1);
                }
                seq_d[full_name] = fa_len;
            } else {
                fa_len = seq_d.at(full_name);
            }
            if (parsed_name == sequence_name && sequence->getSequenceLength() != fa_len) {
                cerr << "[halUnclip]: Sequence " << full_name << " has len=" << fa_len << " in fasta but len=" << sequence->getSequenceLength() << " in hal" << endl;
                exit(1);
            }
            if (parsed_name != sequence_name && sequence->getSequenceLength() > fa_len) {
                cerr << "[halUnclip]: Sequence " << sequence->getFullName() << " has len=" << fa_len << " in fasta but len=" << sequence->getSequenceLength() << " in hal" << endl;
                exit(1);
            }

            frag_map[parsed_name].push_back(sequence);
        }

        // pass 2: compute the dimensions for each base sequence
        for (auto& nf : frag_map) {
            const string& base_name = nf.first;
            string full_name = genome->getName() + "." + base_name;
            vector<const Sequence*>& frags = nf.second;
            size_t fa_len;
            fa_len = seq_d.at(full_name);
            // sort the fragments by start position
            map<size_t, const Sequence*> start_to_frag;
            for (const Sequence* frag : frags) {
                int64_t start = -1;
                string parsed_name = target_set.count(name) ? parse_subpath_name(frag->getName(), &start) : frag->getName();
                if (start == -1) {
                    start = 0;
                    assert(frags.size() == 1);
                }
                start_to_frag[start] = frag;
            }

            // count the top segments
            size_t top = 0;
            // count the gaps (separate counter just for debugging)
            size_t gaps = 0;
            if (start_to_frag.begin()->first > 0) {
                // gap in front
                ++gaps;
            }
            for (auto i = start_to_frag.begin(); i != start_to_frag.end(); ++i) {
                auto next = i;
                ++next;
                if (next != start_to_frag.end()) {
                    if (i->first + i->second->getSequenceLength() < next->first) {
                        // gap in middle
                        ++gaps;
                    }
                }
                top += i->second->getNumTopSegments();                
            }
            if (start_to_frag.rbegin()->first + start_to_frag.rbegin()->second->getSequenceLength() < fa_len) {
                // gap in back
                ++gaps;
            }
            dimensions.push_back(Sequence::Info(base_name, fa_len, top + gaps, 0));
        }
        
        alignment->closeGenome(genome);
    }

    return dim_map;
}

static void copy_and_fill(AlignmentConstPtr in_alignment, AlignmentPtr out_alignment, const unordered_map<string, size_t>& seq_dims,
                          const unordered_set<string>& target_set, bool progress) {
    
    const Genome* in_root_genome = in_alignment->openGenome(in_alignment->getRootName());
    Genome* out_root_genome = out_alignment->openGenome(in_alignment->getRootName());
    
    vector<string> names = in_alignment->getChildNames(in_alignment->getRootName());
    // with a lot of children, the bottom segments are unweildy.  they play havoc with default settings (chunk=1000 is too small)
    // and are terribly slow even with tuning (except inmemory).  so we load up everything we need in this structure in memory
    // so that the bottom segments can be set in a single pass
    vector<vector<hal_index_t>> old_to_new_tsai_vec(names.size());

    for (const string& name : names) {
        if (progress) {
            cerr << "[halUnclip]: Copying segments of " << name << flush;
        }

        const Genome* in_genome = in_alignment->openGenome(name);
        Genome* out_genome = out_alignment->openGenome(name);
        hal_index_t out_child_no = out_root_genome->getChildIndex(out_genome);
        hal_index_t in_child_no = in_root_genome->getChildIndex(in_genome);
        assert(in_child_no == out_child_no);

        // map base name to sequence fragments
        // todo: same thing done in above funciton -- generalize?
        unordered_map<string, vector<const Sequence*>> frag_map;

        // pass 1, map all hal sequences back to their base name and check that they correspond to a fasta sequence
        if (progress) {
            cerr << " [pass 1]" << flush;
        }
        for (SequenceIteratorPtr seqIt = in_genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
            const Sequence *sequence = seqIt->getSequence();
            string sequence_name = sequence->getName();
            string parsed_name = target_set.count(name) ? parse_subpath_name(sequence_name) : sequence_name;
            frag_map[parsed_name].push_back(sequence);
        }

        // pass 2, copy each sequence fragment by fragment
        if (progress) {
            cerr << " [pass 2]" << flush;
        }
        vector<hal_index_t>& old_to_new_tsai = old_to_new_tsai_vec[out_child_no];
        old_to_new_tsai.resize(in_genome->getNumTopSegments(), NULL_INDEX);
    
        for (auto& nf : frag_map) {
            const string& base_name = nf.first;
            vector<const Sequence*>& frags = nf.second;
            
            // sort the fragments by start position
            map<size_t, const Sequence*> start_to_frag;
            for (const Sequence* frag : frags) {
                int64_t start = -1;
                string parsed_name = target_set.count(name) ? parse_subpath_name(frag->getName(), &start) : frag->getName();
                if (start == -1) {
                    start = 0;
                    assert(frags.size() == 1);
                }
                start_to_frag[start] = frag;
            }

            // the one output sequence that corresponds to the list of fragments in the input
            Sequence* out_sequence = out_genome->getSequence(base_name);
            assert(out_sequence != nullptr);
            TopSegmentIteratorPtr out_top = out_sequence->getTopSegmentIterator();
            TopSegment* ts;

            int64_t cur_pos = 0; //position in out_sequence
            int64_t out_start = out_sequence->getStartPosition(); //offset needed when setting coorindatesin out_top

            // visit the ordered input sequence fragments that correspond to out_sequence
            for (auto i = start_to_frag.begin(); i != start_to_frag.end(); ++i) {
                const Sequence* in_sequence_frag = i->second;
                int64_t frag_start = i->first;
                if (frag_start > cur_pos) {
                    // need to add a gap *before* this fragment
                    ts = out_top->tseg();
                    ts->setCoordinates(cur_pos + out_start, frag_start - cur_pos);
                    ts->setParentIndex(NULL_INDEX);
                    ts->setNextParalogyIndex(NULL_INDEX);
                    ts->setBottomParseIndex(NULL_INDEX);
#ifdef debug
                    cerr << "cur_pos=" << cur_pos << flush;
#endif
                    cur_pos += ts->getLength();
#ifdef debug
                    cerr << " after adding start gap cur_pos=" << cur_pos << " (frag name=" << in_sequence_frag->getName() << " fragstart=" << frag_start << ")" << endl;
#endif
                    out_top->toRight();
                }
#ifdef debug
                cerr << "frag " << in_sequence_frag->getFullName() << " has " << in_sequence_frag->getNumTopSegments() << " topsegs which will map to range "
                     << out_sequence->getTopSegmentIterator()->tseg()->getArrayIndex() << " - "
                     << (out_sequence->getTopSegmentIterator()->tseg()->getArrayIndex() + in_sequence_frag->getNumTopSegments()) << endl;
#endif
                // copy the fragment.  note that the ancestor coordinates haven't changed
                // any, so those coordinates can go directly
                TopSegmentIteratorPtr frag_top = in_sequence_frag->getTopSegmentIterator();
                size_t frag_top_count = in_sequence_frag->getNumTopSegments();
                for (size_t frag_top_i = 0; frag_top_i < frag_top_count; ++frag_top_i) {
                    ts = out_top->tseg();
                    ts->setCoordinates(out_start + cur_pos, frag_top->tseg()->getLength());
                    ts->setParentIndex(frag_top->tseg()->getParentIndex());
                    ts->setParentReversed(frag_top->tseg()->getParentReversed());

                    // set the bad value from input alignment, to be update later when we have map
                    ts->setNextParalogyIndex(frag_top->tseg()->getNextParalogyIndex());
                    ts->setBottomParseIndex(NULL_INDEX);
#ifdef debug
                    cerr << "cur_pos=" << cur_pos << flush;
#endif
                    cur_pos += ts->getLength();
#ifdef debug
                    cerr << " after adding frag_ts " << frag_top_i << " cur_pos=" << cur_pos << endl;
#endif
                    old_to_new_tsai[frag_top->tseg()->getArrayIndex()] = ts->getArrayIndex();
                    frag_top->toRight();
                    out_top->toRight();
                }            
            }
            if (cur_pos < (int64_t)out_sequence->getSequenceLength()) {
                // needto add a gap *after* the last fragment
                ts = out_top->tseg();
                ts->setCoordinates(out_start + cur_pos, (int64_t)out_sequence->getSequenceLength() - cur_pos);
                ts->setParentIndex(NULL_INDEX);
                ts->setNextParalogyIndex(NULL_INDEX);
                ts->setBottomParseIndex(NULL_INDEX);
#ifdef debug
                cerr << "cur_pos="<< cur_pos << flush;
#endif
                cur_pos += ts->getLength();
#ifdef debug
                cerr << " after adding end gap cur_pos=" << cur_pos << endl;
#endif
                out_top->toRight();
            }
            if (cur_pos != (int64_t)out_sequence->getSequenceLength()) {
                cerr << "[halUnclip]: sanity check fail for sequence " << name << "." << base_name << ".  The offset after conversion is "
                     << cur_pos << " which is different than the sequence length of " << out_sequence->getSequenceLength() << endl
                     << "[halUnclip]: the fragments are\n";
                for (size_t i = 0; i < frags.size(); ++i) {
                    const Sequence* in_sequence_frag = frags[i];
                    cerr << "     " << in_sequence_frag->getName() << " len=" << in_sequence_frag->getSequenceLength() << endl;                    
                }
            }
            assert(cur_pos == (int64_t)out_sequence->getSequenceLength());            
            assert(out_top->getArrayIndex() == out_sequence->getTopSegmentIterator()->getArrayIndex() + (int64_t)out_sequence->getNumTopSegments());
        }

        //pass 3: set the paralogy indexes
        if (progress) {
            cerr << " [pass 3]" << endl;
        }        
        TopSegment* ts;
        for (TopSegmentIteratorPtr out_topit = out_genome->getTopSegmentIterator(); !out_topit->atEnd(); out_topit->toRight()) {
            ts = out_topit->tseg();
            if (ts->hasNextParalogy()) {
                ts->setNextParalogyIndex(old_to_new_tsai[ts->getNextParalogyIndex()]);
            }
        }

        in_alignment->closeGenome(in_genome);
        out_alignment->closeGenome(out_genome);
    }

    // copy the root
    if (progress) {
        cerr << "[halUnclip]: Copying root segments" << endl;
    }
    BottomSegmentIteratorPtr in_botit = in_root_genome->getBottomSegmentIterator();
    BottomSegmentIteratorPtr out_botit = out_root_genome->getBottomSegmentIterator();
    assert(in_root_genome->getNumBottomSegments() == out_root_genome->getNumBottomSegments());
    assert(in_root_genome->getNumChildren() == out_root_genome->getNumChildren());
    size_t num_bottom = in_root_genome->getNumBottomSegments();
    size_t num_children = in_root_genome->getNumChildren();
    for (size_t i = 0; i < num_bottom; ++i) {
        BottomSegment* in_bs = in_botit->bseg();
        BottomSegment* out_bs = out_botit->bseg();
        out_bs->setCoordinates(in_bs->getStartPosition(), in_bs->getLength());
        for (size_t j = 0; j < num_children; ++j) {
            // everything's the same except the child index, which gets mapped via old_to_new_tsai_vec
            hal_index_t in_ci = in_bs->getChildIndex(j);
            hal_index_t out_ci = in_ci != NULL_INDEX ? old_to_new_tsai_vec[j][in_ci] : in_ci;
            out_bs->setChildIndex(j, out_ci);
            out_bs->setChildReversed(j, in_bs->getChildReversed(j));
        }
        out_bs->setTopParseIndex(NULL_INDEX);
        in_botit->toRight();
        out_botit->toRight();
    }
}

// go in and rewerite the sequences from the fasta
void add_fasta_sequences(AlignmentConstPtr in_alignment, AlignmentPtr out_alignment, const string& seqfile_path, const unordered_set<string>& target_set, bool progress) {
    ifstream seqfile(seqfile_path);
    if (!seqfile) {
        cerr << "[halUnclip]: Unable to open seqfile: " << seqfile_path << endl;
        exit(1);
    }

    string buffer;
    set<string> done_set;
    while (getline(seqfile, buffer)) {
        vector<string> toks = split_delims(buffer, " \t");
        if (toks.size() == 2) {
            string name = toks[0];
            string fa_path = toks[1];
            if (target_set.count(name)) {
                done_set.insert(name);
                if (progress) {
                    cerr << "[halUnclip]: Loading fasta for " << name << " ... " << flush;
                }
                unordered_map<string, pair<size_t, string>> fa_info = read_fasta(fa_path);
                if (progress) {
                    cerr << "and setting dna strings in output genome" << endl;
                }
                Genome* genome = out_alignment->openGenome(name);
                assert(genome != nullptr);
                for (auto& fi : fa_info) {
                    Sequence* sequence = genome->getSequence(fi.first);
                    if (sequence != nullptr) {
                        assert(sequence->getSequenceLength() == fi.second.first);
                        sequence->setString(fi.second.second);
                    }
                }
            }
        }
    }

    // if there's no _sub sequences found, a genome is allowed to not be in the sequence map
    // this is generally the case for the root, but could be the minigraph contigs
    vector<string> names = in_alignment->getChildNames(in_alignment->getRootName());
    names.push_back(in_alignment->getRootName());
    for (const string& name : names) {
        if (!done_set.count(name)) {
            if (progress) {
                cerr << "[halUnclip]: Directly copying dna strings for " << name << endl;
            };
            const Genome* in_genome = in_alignment->openGenome(name);
            Genome* out_genome = out_alignment->openGenome(name);
            for (SequenceIteratorPtr seqIt = in_genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
                const Sequence* in_sequence = seqIt->getSequence();
                Sequence* out_sequence = out_genome->getSequence(in_sequence->getName());
                in_sequence->getString(buffer);
                out_sequence->setString(buffer);
            }
            if (name != in_alignment->getRootName()) {
                in_alignment->closeGenome(in_genome);
                out_alignment->closeGenome(out_genome);
            }
        }
    }
}


// root->leaf alignments are consistent
static void validate_alignments(AlignmentConstPtr in_alignment, AlignmentPtr out_alignment) {

    validateAlignment(out_alignment.get());

    const Genome* in_root_genome = in_alignment->openGenome(in_alignment->getRootName());
    Genome* out_root_genome = out_alignment->openGenome(in_alignment->getRootName());
    assert(in_root_genome->getNumBottomSegments() == out_root_genome->getNumBottomSegments());
    assert(in_root_genome->getNumChildren() == out_root_genome->getNumChildren());
    // we go by genome (instead of segment) to hopefully be cache-friendlier
    for (size_t j = 0; j < in_root_genome->getNumChildren(); ++j) {
        const Genome* in_genome = in_root_genome->getChild(j);
        Genome* out_genome = out_root_genome->getChild(j);
        BottomSegmentIteratorPtr in_botit = in_root_genome->getBottomSegmentIterator();
        BottomSegmentIteratorPtr out_botit = out_root_genome->getBottomSegmentIterator();
        TopSegmentIteratorPtr in_topit = in_genome->getTopSegmentIterator();
        TopSegmentIteratorPtr out_topit = out_genome->getTopSegmentIterator();
        for (size_t i = 0; i < in_genome->getNumBottomSegments(); ++i) {
            in_topit->toChild(in_botit, j);
            out_topit->toChild(out_botit, j);
            
            string s1, s2;
            if (j == 0) {
                in_botit->getString(s1);
                out_botit->getString(s2);
                assert(s1 == s2);
            }
            in_topit->getString(s1);
            out_topit->getString(s2);
            assert(s1 == s2);

            string in_seq_name = in_topit->tseg()->getSequence()->getName();
            string out_seq_name = out_topit->tseg()->getSequence()->getName();
            int64_t start;
            string in_base_name = parse_subpath_name(in_seq_name, &start);
            assert(in_base_name == out_seq_name);
            assert(in_topit->getReversed() == out_topit->getReversed());
            if (!in_topit->getReversed()) {
                // punt on reverse check for now
                assert(in_topit->getStartPosition() + start == out_topit->getStartPosition());
            }

            in_botit->toRight();
            out_botit->toRight();
        }
        in_alignment->closeGenome(in_genome);
        out_alignment->closeGenome(out_genome);        
    }
}

int main(int argc, char** argv) {
    CLParser optionsParser(WRITE_ACCESS);
    initParser(&optionsParser);
    string in_hal_path;
    string out_hal_path;
    string seqfile_path;
    string target_genomes;
    bool progress;
    bool validate;
    try {
        optionsParser.parseOptions(argc, argv);
        in_hal_path = optionsParser.getArgument<string>("inFile");
        seqfile_path = optionsParser.getArgument<string>("seqFile");
        out_hal_path = optionsParser.getArgument<string>("outFile");
        target_genomes = optionsParser.getOption<string>("targetGenomes");
        progress = optionsParser.getFlag("progress");
        validate = optionsParser.getFlag("validate");
    }
    catch(exception& e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }

    // load the input genome
    if (progress) {
        cerr << "[halUnclip]: Opening input alignment" << endl;
    }
    AlignmentConstPtr in_alignment(openHalAlignment(in_hal_path, &optionsParser, READ_ACCESS));

    // and the output genome
    if (progress) {
        cerr << "[halUnclip]: Creating output alignment object" << endl;
    }    
    AlignmentPtr out_alignment(openHalAlignment(out_hal_path, &optionsParser, READ_ACCESS | WRITE_ACCESS | CREATE_ACCESS));

    // check the targets, defaulting to all leaves
    vector<string> target_names;
    if (target_genomes != "\"\"") {
        target_names = chopString(target_genomes, ",");
        for (const string& name : target_names) {
            const Genome* genome = in_alignment->openGenome(name);
            if (genome == nullptr) {
                cerr << "[halUnclip]: Target genome " << name << " not present in input HAL" << endl;
                exit(1);
            }
            in_alignment->closeGenome(genome);
        }
    } else {
        target_names = in_alignment->getChildNames(in_alignment->getRootName());
    }
    unordered_set<string> target_set(target_names.begin(), target_names.end());

    // and load the fasta sequence sizes from the seqfile
    if (progress) {
        cerr << "[halUnclip]: Reading fasta dimensions from seqfile" << endl;
    }
    unordered_map<string, size_t> seq_dims = get_dimensions_from_seqfile(seqfile_path, target_set);

    if (progress) {
        cerr << "[halUnclip]: Computing new hal dimensions" << endl;
    }
    unordered_map<string, vector<Sequence::Info>> dimensions = get_filled_dimensions(in_alignment, seq_dims, target_set, progress);

    // set up the size of each genome, staring with the root
    string root_name = in_alignment->getRootName();
    Genome* root_genome = out_alignment->addRootGenome(root_name);
    // important to visit these in order, so child indexes are presesrved
    vector<string> leaf_names = in_alignment->getChildNames(root_name);
    for (const string& leaf_name : leaf_names) {
        vector<Sequence::Info>& leaf_dims = dimensions.at(leaf_name);
        Genome* leaf_genome = out_alignment->addLeafGenome(leaf_name, root_name, 1);
        leaf_genome->setDimensions(leaf_dims);
        if (progress) {
            cerr << "[halUnclip]: Adding leaf genome " << leaf_name << " with length " << leaf_genome->getSequenceLength() << " and " << leaf_genome->getNumTopSegments() << " top segments" << endl;
        }
    }

    // important to set root dimensions after adding leaves so bottom segments have right number of slots
    root_genome->setDimensions(dimensions.at(root_name));
    if (progress) {
        cerr << "[halUnclip]: Adding root genome " << root_name << " with length " << root_genome->getSequenceLength() << " and " << root_genome->getNumBottomSegments() << " bottom segments" << endl;
    }    
    
    // copy over the filled graph
    if (progress) {
        cerr << "[halUnclip]: Copying and filling the graph" << endl;
    }
    copy_and_fill(in_alignment, out_alignment, seq_dims, target_set, progress);

    // add back the fasta sequences
    if (progress) {
        cerr << "[halUnclip]: Adding fasta sequences" << endl;
    }    
    add_fasta_sequences(in_alignment, out_alignment, seqfile_path, target_set, progress);

    if (validate) {
        if (progress) {
            cerr << "[halUnclip]: Validating alignment" << endl;
        }
        validate_alignments(in_alignment, out_alignment);
    }    
            
    if (progress) {
        cerr << "[halUnclip]: Writing output alignment" << endl;
    }
     
    return 0;
}

