/**
 * @file context_positions.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implements a class to collect context positions
 * @version 0.6
 * @date 2023-07-29
 * 
 * @copyright Copyright (c) 2023
 * 
 * MIT License
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __RACES_CONTEXT_POSITIONS__
#define __RACES_CONTEXT_POSITIONS__

#include <map>
#include <set>
#include <memory>
#include <vector>
#include <fstream>
#include <type_traits>

#include "archive.hpp"
#include "fasta_utils.hpp"      // IO::FASTA::EnsemblSeqNameDecoder
#include "fragment.hpp"         // Passengers::GenomicRegion
#include "basic_IO.hpp"         // IO::get_stream_size
#include "context.hpp"          // Passengers::ExtendedContextAutomaton

#include "progress_bar.hpp"


namespace Races
{

namespace Passengers
{

template<typename ABSOLUTE_GENOMIC_POSITION = uint32_t,
         std::enable_if_t<std::is_integral_v<ABSOLUTE_GENOMIC_POSITION>
                          && std::is_unsigned_v<ABSOLUTE_GENOMIC_POSITION>, bool> = true>
struct ContextPositions
{
    using AbsoluteGenomicPosition = ABSOLUTE_GENOMIC_POSITION;
protected:
    using ContextPositionMap = std::map<MutationalContext, std::vector<AbsoluteGenomicPosition> >;

    std::shared_ptr<ContextPositionMap> context2pos;                  //!< the context-genomic positions map
    std::map<AbsoluteGenomicPosition, ChromosomeId> abs_pos2chr;      //!< the absolute genomic position-chromosome id map

    AbsoluteGenomicPosition genome_size;  //!< the genome size

    /**
     * @brief Find the mutational contexts of a FASTA sequence and save their positions in the map
     * 
     * This method finds the mutational contexts of a FASTA sequence and saves their absolute positions
     * in `ContextPositions::context2pos`.
     * 
     * @param[in,out] fasta_stream is the FASTA stream pointing at the first nucleotide of the considered sequence
     * @param[in,out] progress_bar is the progress bar
     */
    void find_contexts_in_seq(std::ifstream& fasta_stream, UI::ProgressBar* progress_bar)
    {
        auto filesize = IO::get_stream_size(fasta_stream);
        if (progress_bar != nullptr) {
            progress_bar->set_progress(static_cast<uint8_t>(100*fasta_stream.tellg()/filesize));
        }

        ExtendedContextAutomaton c_automata;

        char last_char{'A'};
        while (last_char != '>' && !fasta_stream.eof()) {
            fasta_stream.get(last_char);

            if (c_automata.update_state(last_char)) {
                if (c_automata.read_a_context()) {
                    (*context2pos)[c_automata.get_context()].emplace_back(genome_size);
                }
                ++genome_size;
            }

            // update progress bar once every 2^22-1 nucleotides
            if ((genome_size&0x3FFFFF)==0 && progress_bar != nullptr) {
                progress_bar->set_progress(static_cast<uint8_t>(100*fasta_stream.tellg()/filesize));
            }
        }
    }

    /**
     * @brief Skip nucleotides in a FASTA stream up to a specified position in the current sequence
     * 
     * @param[in,out] fasta_stream is the input FASTA stream
     * @param[in,out] c_automata is an extended context automaton
     * @param[in,out] current_pos is the current position
     * @param[in] aimed_position is the aimed position
     */
    void skip_to(std::ifstream& fasta_stream, ExtendedContextAutomaton& c_automata, 
                 ChrPosition& position, ChrPosition aimed_position)
    {
        char last_char{'N'};
        while (last_char != '>' && !fasta_stream.eof() && position < aimed_position) {
            fasta_stream.get(last_char);

            if (c_automata.update_state(last_char)) {
                ++position;
                ++genome_size;
            }
        }
    }

    /**
     * @brief Skip all the remaining nucleotides of the current sequence in a FASTA stream
     * 
     * This method discharge the remaining part of the  current sequence, but it does consider
     * it as part of the genome.
     * 
     * @param[in,out] fasta_stream is the input FASTA stream
     */
    void skip_to_next_seq(std::ifstream& fasta_stream)
    {
        char in_char{'N'};
        while (in_char != '>' && !fasta_stream.eof())
        {
            fasta_stream.get(in_char);

            if (ExtendedContextAutomaton::is_a_nucleotide(in_char)) {
                ++genome_size;
            }
        }
    }

    /**
     * @brief Skip all the remaining nucleotides of the current sequence in a FASTA stream
     * 
     * This method discharge the current sequence and does not consider it as part of the
     * genome.
     * 
     * @param[in,out] fasta_stream is the input FASTA stream
     */
    static void discharge_seq(std::ifstream& fasta_stream)
    {
        char in_char{'N'};
        while (in_char != '>' && !fasta_stream.eof())
        {
            fasta_stream.get(in_char);
        }
    }

    /**
     * @brief Find the mutational contexts in parts of a FASTA sequence and save their positions in the map
     * 
     * This method finds the mutational contexts laying in the specified regions of a FASTA sequence and 
     * saves their absolute positions in `ContextPositions::context2pos`.
     * 
     * @param[in,out] fasta_stream is the FASTA stream pointing at the first nucleotide of the considered sequence
     * @param[in] genomic_regions is the set of regions to be considered when searching mutational contexts
     * @param[in,out] progress_bar is the progress bar
     */
    void find_contexts_in_seq(std::ifstream& fasta_stream, const std::set<GenomicRegion>& genomic_regions,
                              UI::ProgressBar* progress_bar)
    {
        auto filesize = Races::IO::get_stream_size(fasta_stream);
        if (progress_bar != nullptr) {
            progress_bar->set_progress(static_cast<uint8_t>(100*fasta_stream.tellg()/filesize));
        }

        ExtendedContextAutomaton c_automata;

        char last_char{'N'};
        ChrPosition chr_pos{0};
        auto region_it = genomic_regions.begin();
        while (last_char != '>' && !fasta_stream.eof() && region_it != genomic_regions.end()) {
            skip_to(fasta_stream, c_automata, chr_pos, region_it->get_initial_position()-1);

            while (last_char != '>' && !fasta_stream.eof() && chr_pos!=region_it->get_final_position()+1) {
                fasta_stream.get(last_char);

                if (c_automata.update_state(last_char)) {
                    if (c_automata.read_a_context()) {
                        (*context2pos)[c_automata.get_context()].emplace_back(genome_size);
                    }
                    ++chr_pos;
                    ++genome_size;
                }

                // update progress bar once every 2^22-1 nucleotides
                if ((genome_size&0x3FFFFF)==0 && progress_bar != nullptr) {
                    progress_bar->set_progress(static_cast<uint8_t>(100*fasta_stream.tellg()/filesize));
                }
            }

            ++region_it;
        }

        if (last_char != '>' && !fasta_stream.eof()) {
            skip_to_next_seq(fasta_stream);
        }
    }

    /**
     * @brief Split a set of genomic regions by chromosome id
     * 
     * @param[in] genomic_regions is the set of genomic region to be split 
     * @return a map that associates a chromosome id to the the set of genomic regions
     *     laying in the corresponding chromosome
     */
    static std::map<ChromosomeId, std::set<GenomicRegion> > split_by_chromosome_id(const std::set<GenomicRegion>* genomic_regions)
    {
        std::map<ChromosomeId, std::set<GenomicRegion> > split;

        if (genomic_regions != nullptr) {
            for (const auto& genomic_region: *genomic_regions) {
                split[genomic_region.get_chromosome_id()].insert(genomic_region);
            }
        }
        return split;
    }

    /**
     * @brief Reset the object and find the mutational contexts in a FASTA stream
     * 
     * This method resets the current object, finds the mutational contexts in the chromosome sequences
     * contained in a FASTA stream, and saves their absolute positions in the object member 
     * `ContextPositions::context2pos`. The chromosome sequences are identified by their names according 
     * to the method template parameter.
     * At the same time, this method stores the absolute positions of each chromosome in the member 
     * `ContextPositions::abs_pos2chromosome`.
     * 
     * @tparam SEQUENCE_NAME_DECODER is a sequence name decoder inherited from `Races::IO::FASTA::SeqNameDecoder` 
     * @param[in,out] fasta_stream is the genome FASTA stream
     * @param[in] genomic_regions is the vector of the regions to be processed
     * @param[in,out] progress_bar is the progress bar
     */
    template<typename SEQUENCE_NAME_DECODER = Races::IO::FASTA::EnsemblSeqNameDecoder,
             std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SeqNameDecoder, SEQUENCE_NAME_DECODER>, bool> = true >
    void reset_with(std::ifstream& fasta_stream, const std::set<GenomicRegion>* genomic_regions=nullptr,
                    UI::ProgressBar* progress_bar=nullptr)
    {
        if (!fasta_stream.good()) {
            throw std::runtime_error("the stream is not readable");
        }

        context2pos = std::make_shared<ContextPositionMap>();
        abs_pos2chr.clear();
        genome_size = 0;

        auto regions_by_chr = split_by_chromosome_id(genomic_regions);

        SEQUENCE_NAME_DECODER seq_decoder;
        skip_to_next_seq(fasta_stream);

        while (!fasta_stream.eof()) {
            std::string sequence_title;
            getline(fasta_stream, sequence_title);

            GenomicRegion chr_region;

            // if the sequence is a chromosome
            if (seq_decoder.get_chromosome_region(sequence_title, chr_region)) {

                if (progress_bar != nullptr) {
                    progress_bar->set_message("Processing chr. " + GenomicPosition::chrtos(chr_region.get_chromosome_id()));
                }
                abs_pos2chr[genome_size+1] = chr_region.get_chromosome_id();

                if (genomic_regions == nullptr) {
                    find_contexts_in_seq(fasta_stream, progress_bar);
                } else {
                    auto regions_in = regions_by_chr.find(chr_region.get_chromosome_id());

                    // if one of the selected region belongs to the current chromosome
                    if (regions_in != regions_by_chr.end()) {
                        find_contexts_in_seq(fasta_stream, regions_in->second, progress_bar);
                    } else {
                        skip_to_next_seq(fasta_stream);
                    }
                }
            } else {
                discharge_seq(fasta_stream);
            }
        }

        if (progress_bar != nullptr) {
            progress_bar->set_progress(100, "Context index build");
        }
    }
public:

    /**
     * @brief The empty constructor
     */
    ContextPositions():
        context2pos(std::make_shared<ContextPositionMap>()), genome_size(0)
    {}

    /**
     * @brief Find the context positions in a FASTA stream
     * 
     * @tparam SEQUENCE_NAME_DECODER is a sequence name decoder inherited from `Races::IO::FASTA::SeqNameDecoder`
     * @param[in,out] genome_fasta_stream is a input file stream in FASTA format
     * @param[in,out] progress_bar is the progress bar 
     * @return the positions of the contexts in the FASTA sequences whose name 
     *      corresponds to a chromosome according to `SEQUENCE_NAME_DECODER`
     */
    template<typename SEQUENCE_NAME_DECODER = Races::IO::FASTA::EnsemblSeqNameDecoder,
             std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SeqNameDecoder, SEQUENCE_NAME_DECODER>, bool> = true >
    static ContextPositions find_contexts_in(std::ifstream& genome_fasta_stream, UI::ProgressBar* progress_bar=nullptr)
    {
        ContextPositions c_positions;

        c_positions.reset_with<SEQUENCE_NAME_DECODER>(genome_fasta_stream, nullptr, progress_bar);

        return c_positions;
    }

    /**
     * @brief Find the context positions in a FASTA file
     * 
     * @tparam SEQUENCE_NAME_DECODER is a sequence name decoder inherited from `Races::IO::FASTA::SeqNameDecoder`
     * @param[in] genome_fasta genome_fasta is a path of a FASTA file
     * @param[in,out] progress_bar is the progress bar 
     * @return the positions of the contexts in the FASTA file whose name 
     *      corresponds to a chromosome according to `SEQUENCE_NAME_DECODER`
     */
    template<typename SEQUENCE_NAME_DECODER = Races::IO::FASTA::EnsemblSeqNameDecoder,
             std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SeqNameDecoder, SEQUENCE_NAME_DECODER>, bool> = true >
    static ContextPositions find_contexts_in(const std::filesystem::path& genome_fasta, UI::ProgressBar* progress_bar=nullptr)
    {
        std::ifstream genome_fasta_stream(genome_fasta);

        if (!genome_fasta_stream.good()) {
            std::ostringstream oss;
            
            oss << "\"" << genome_fasta << "\" does not exist"; 
            throw std::runtime_error(oss.str());
        }        
        
        return find_contexts_in<SEQUENCE_NAME_DECODER>(genome_fasta_stream, progress_bar);
    }

    /**
     * @brief Find the context positions in some genomic fragments of a FASTA stream
     * 
     * @tparam SEQUENCE_NAME_DECODER is a sequence name decoder inherited from `Races::IO::FASTA::SeqNameDecoder`
     * @param[in,out] genome_fasta_stream is a input file stream in FASTA format
     * @param[in] genomic_regions is the set of the genomic regions in which contexts will be searched 
     * @param[in,out] progress_bar is the progress bar 
     * @return the positions of the contexts in the regions `genomic_regions` of the FASTA 
     *      sequences whose name corresponds to a chromosome according to `SEQUENCE_NAME_DECODER`
     */
    template<typename SEQUENCE_NAME_DECODER = Races::IO::FASTA::EnsemblSeqNameDecoder,
             std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SeqNameDecoder, SEQUENCE_NAME_DECODER>, bool> = true >
    static ContextPositions find_contexts_in(std::ifstream& genome_fasta_stream, const std::set<GenomicRegion>& genomic_regions,
                                             UI::ProgressBar* progress_bar=nullptr)
    {
        ContextPositions c_positions;

        c_positions.reset_with<SEQUENCE_NAME_DECODER>(genome_fasta_stream, &genomic_regions, progress_bar);

        return c_positions;
    }

    /**
     * @brief Find the context positions in some genomic fragments of a FASTA file
     * 
     * @tparam SEQUENCE_NAME_DECODER is a sequence name decoder inherited from `Races::IO::FASTA::SeqNameDecoder`
     * @param[in] genome_fasta genome_fasta is a path of a FASTA file
     * @param[in] genomic_regions is the set of the genomic regions in which contexts will be searched 
     * @param[in,out] progress_bar is the progress bar 
     * @return the positions of the contexts in the regions `genomic_regions` of the FASTA 
     *      sequences whose name corresponds to a chromosome according to `SEQUENCE_NAME_DECODER`
     */
    template<typename SEQUENCE_NAME_DECODER = Races::IO::FASTA::EnsemblSeqNameDecoder,
             std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SeqNameDecoder, SEQUENCE_NAME_DECODER>, bool> = true >
    static ContextPositions find_contexts_in(const std::filesystem::path& genome_fasta, const std::set<GenomicRegion>& genomic_regions,
                                             UI::ProgressBar* progress_bar=nullptr)
    {
        std::ifstream genome_fasta_stream(genome_fasta);

        if (!genome_fasta_stream.good()) {
            std::ostringstream oss;
            
            oss << "\"" << genome_fasta << "\" does not exist"; 
            throw std::runtime_error(oss.str());
        }      

        return find_contexts_in<SEQUENCE_NAME_DECODER>(genome_fasta_stream, genomic_regions, progress_bar);
    }

    /**
     * @brief Get the context positions
     * 
     * @return a constant reference to simulator context positions
     */
    inline const std::map<MutationalContext, std::vector<ABSOLUTE_GENOMIC_POSITION> >& get_context_positions() const
    {
        return *context2pos;
    }

    /**
     * @brief Get the absolute genomic positions of a context
     * 
     * @param context is the searched context
     * @return the vector of the absolute genomic positions 
     */
    inline const std::vector<ABSOLUTE_GENOMIC_POSITION>& operator[](const MutationalContext& context) const
    {
        return context2pos->at(context);
    }

    /**
     * @brief Get the regions of the chromosomes
     * 
     * @return the regions of the chromosomes
     */
    std::vector<GenomicRegion> get_chromosome_regions() const
    {
        std::vector<GenomicRegion> chr_regions;

        if (abs_pos2chr.size()>0) {
            auto it = abs_pos2chr.begin();
            auto next = it;

            while (++next != abs_pos2chr.end()) {
                chr_regions.emplace_back(it->second, static_cast<GenomicRegion::Length>(next->first-it->first));
                it = next;
            }

            chr_regions.emplace_back(it->second, static_cast<GenomicRegion::Length>(genome_size-it->first+1));
        }

        return chr_regions;
    }

    /**
     * @brief Turn an absolute genomic position into the corresponding genomic position
     * 
     * @param abs_position is an absolute genomic position
     * @return the genomic position corresponding to `abs_position`  
     */
    GenomicPosition get_genomic_position(const ABSOLUTE_GENOMIC_POSITION& abs_position) const
    {
        auto it = abs_pos2chr.upper_bound(abs_position);
        --it;

        return GenomicPosition(it->second, static_cast<ChrPosition>(abs_position-it->first));
    }

    /**
     * @brief Get the genome size
     * 
     * @return get the genome size 
     */
    inline const ABSOLUTE_GENOMIC_POSITION& get_genome_size() const
    {
        return genome_size;
    }

    /**
     * @brief Save a simulator in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    inline void save(ARCHIVE& archive) const
    {
        uint8_t abs_pos_size = sizeof(ABSOLUTE_GENOMIC_POSITION);

        archive & abs_pos_size
                & *context2pos
                & abs_pos2chr
                & genome_size;
    }

    /**
     * @brief Load a simulator from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded simulator
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static ContextPositions load(ARCHIVE& archive)
    {
        ContextPositions c_positions;
        uint8_t abs_pos_size;

        archive & abs_pos_size;

        if (abs_pos_size != sizeof(ABSOLUTE_GENOMIC_POSITION)) {
            std::ostringstream oss;

            oss << "Absolute position size (" << sizeof(ABSOLUTE_GENOMIC_POSITION) 
                << " bytes) does not match file object one (" << static_cast<size_t>(abs_pos_size) 
                << " bytes)." ;

            throw std::runtime_error(oss.str());
        }

        c_positions.context2pos = std::make_shared<ContextPositionMap>();

        archive & *(c_positions.context2pos)
                & c_positions.abs_pos2chr
                & c_positions.genome_size;

        return c_positions;
    }
};

}   // Races

}   // Passengers


#endif // __RACES_CONTEXT_POSITIONS__