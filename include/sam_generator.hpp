/**
 * @file sam_generator.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Defines classes to generate genomic mutations SAM file
 * @version 0.2
 * @date 2023-09-07
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

#ifndef __RACES_SAM_GENERATOR__
#define __RACES_SAM_GENERATOR__

#include <list>
#include <map>
#include <string>
#include <sstream>
#include <filesystem>

#include "genomic_position.hpp"
#include "genome_mutations.hpp"

#include "fasta_reader.hpp"
#include "fasta_utils.hpp"

#include "progress_bar.hpp"

namespace Races
{

namespace Passengers
{

/**
 * @brief The Input/Output namespace for classes related to passenger mutations
 */
namespace IO
{

/**
 * @brief A genomic mutations SAM producer 
 * 
 * @tparam RANDOM_GENERATOR is the random generator type used to sample genome
 */
template<typename RANDOM_GENERATOR=std::mt19937_64>
class SAMGenerator
{
    RANDOM_GENERATOR random_generator;  //!< The random generator

    std::filesystem::path SAM_directory_name;   //!< The SAM output directory
    std::filesystem::path ref_genome_filename;  //!< The reference genome FASTA filename

    size_t read_size;       //!< The produced read size
    size_t num_of_reads;    //!< Number of already placed reads

    /**
     * @brief Chromosome data
     * 
     * This template represents chromosome data. The objects of this class store 
     * chromosome information, i.e., name, header, size, and chromosome identifier,
     * and, depending on the parameter, which must be in the hierarchy of the class
     * `Races::IO::FASTA::SequenceInfo`, may also maintain the chromosome nucleic
     * sequence.
     * 
     * @tparam DATA_TYPE is the base type of the template. It must be a inherited
     *      from `Races::IO::FASTA::SequenceInfo`
     */
    template<typename DATA_TYPE, std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
    struct ChromosomeData : public DATA_TYPE
    {
        ChromosomeId chr_id;    //!< the chromosome id
    };

    struct FilterNonChromosomeSequence : public Races::IO::FASTA::SequenceFilter
    {
        ChromosomeId last_chr_id;   //!< the identifier of the last processed chromosome header

        inline bool operator()(const std::string& header)
        {
            return !Races::IO::FASTA::is_chromosome_header(header, last_chr_id);
        } 
    };

    /**
     * @brief Get the overall allelic size of a list of genomic mutations
     * 
     * This method computes the sum of the allelic sizes of the genomic 
     * mutations in a list.
     * 
     * @param mutations is the list of genomic mutations
     * @param progress_bar is a progress bar
     * @return the sum of the allelic sizes of the genomic mutations in `mutations`
     */
    static size_t get_overall_allelic_size(const std::list<GenomeMutations>& mutations,
                                           UI::ProgressBar& progress_bar)
    {
        progress_bar.set_message("Evaluating allelic size");
        size_t overall_allelic_size{0};
        for (const auto& mutations: mutations) {
            overall_allelic_size += mutations.allelic_size();
            progress_bar.update_elapsed_time();
        }

        return overall_allelic_size;
    }

    /**
     * @brief Check whether all the genomic mutations in a list represent the same genome
     * 
     * This method checks whether all the genomic mutations in a list represent genomes having 
     * the same number of chromosomes and the chromosome sizes are consistent.
     * 
     * @param mutations is a list of genome mutations
     * @param progress_bar is a progress bar
     * @return `true` if and only if all the genomic mutations in `mutations` have the same 
     *          number of chromosomes and, for all pairs of elements in `mutations`, \f$m_1\f$
     *          and \f$m_2\f$, the $i$-th chromosomes in \f$m_1\f$ and \f$m_2\f$ are equal in
     *          size
     */
    static bool represents_same_genome(const std::list<GenomeMutations>& mutations,
                                       UI::ProgressBar& progress_bar)
    {
        if (mutations.size()==0) {
            return true;
        }

        // get the first genome mutations in the list
        const auto& first_mutations = mutations.front();
        
        progress_bar.set_message("Checking input");
        // for all the genome mutations in the list
        for (const auto& cell_mutations : mutations) {

            // if the first element and the current one differ in the number of chromosome
            if (cell_mutations.get_chromosomes().size()!=first_mutations.get_chromosomes().size()) {
                return false;
            }

            // for all the chromosomes in the two genome mutations
            auto cell_it = cell_mutations.get_chromosomes().begin();
            auto first_it = first_mutations.get_chromosomes().begin();
            for (;cell_it != cell_mutations.get_chromosomes().end(); ++cell_it, ++first_it) {
                // if they differ because of their chromosome identifier or because of their 
                // sizes
                if ((cell_it->first != first_it->first)
                        ||(cell_it->second.size()!=first_it->second.size())) {
                    return false;
                }
            }
            progress_bar.update_elapsed_time();
        }

        // otherwise, return true
        return true;
    }

    /**
     * @brief Get the CIGAR string of a genomic region potentially affected by SNVs
     * 
     * The CIGAR string a code that represents the alignment of a sequence over a reference
     * in SAM files (see [1]). This method considers a set of SNVs, a genomic position, 
     * and a size, and it produces the CIGAR code corresponding to a simulated read whose 
     * sequence correspond to that of the reference genome from the specified position with
     * the exception of positions in which the given SNVs were applied.
     * 
     * [1] "Sequence Alignment/Map Format Specification", The SAM/BAM Format Specification 
     *     Working Group, 22 August 2022, http://samtools.github.io/hts-specs/SAMv1.pdf
     * 
     * @param SNVs is a map from genomic positions to the corresponding SNVs
     * @param genomic_position is the initial position of the simulated read
     * @param read_size is the size of the simulated read
     * @return the CIGAR code corresponding to a simulated read having size `read_size` 
     *      whose sequence correspond to that of the sequence in the reference genome 
     *      from `genomic_position` with the exception of the genomic position in which 
     *      the SNVs in `SNVs` were placed
     */
    static std::string get_CIGAR(const std::map<GenomicPosition, SNV>& SNVs,
                                 const GenomicPosition& genomic_position, const size_t& read_size)
    {
        using namespace Races::Passengers;

        const ChrPosition last_position = genomic_position.position+read_size-1;

        std::ostringstream oss;

        // the first position in the last matching sequence
        size_t first_match = 0;

        // for all the SNVs placed on the simulated read
        for (auto snv_it = SNVs.lower_bound(genomic_position);
                snv_it != SNVs.end() && snv_it->second.position<=last_position; ++snv_it) {

            // get the SNV position in the simulated read
            const size_t SNV_pos = snv_it->second.position-genomic_position.position;

            // if we have matched some nucleotides
            if (SNV_pos+1>first_match) {
                // report the number of matched nucleotides
                oss << (SNV_pos-first_match) << "M";
            }

            // report the SNV mismatch
            oss << "1X";

            // the next matching sequence will start from the next 
            // position in the read
            first_match = SNV_pos+1;
        }

        // if the last position was a match
        if (genomic_position.position+first_match <= last_position+1) {
            // report the number of matched nucleotides
            oss << (last_position-genomic_position.position+1-first_match) << "M";
        }

        return oss.str();
    }

    /**
     * @brief Apply a set of SNVs to a fragment of the reference genome
     * 
     * @param nucleic_sequence is a fragment of the reference genome
     * @param SNVs is a map from genomic positions to the corresponding SNVs
     * @param genomic_position is the initial position of `nucleic_sequence` in
     *          the reference genome
     * @return a reference to the modified version of `nucleic_sequence` where 
     *          the SNVs that were placed in the reference genome region 
     *          corresponding to `nucleic_sequence` are applied
     */
    static std::string& apply_SNVs(std::string& nucleic_sequence,
                                   const std::map<GenomicPosition, SNV>& SNVs,
                                   const GenomicPosition& genomic_position)
    {
        using namespace Races::Passengers;

        const ChrPosition last_position = genomic_position.position+nucleic_sequence.size()-1;

        for (auto snv_it = SNVs.lower_bound(genomic_position);
                snv_it != SNVs.end() && snv_it->second.position<=last_position; ++snv_it) {
            const size_t SNV_pos = snv_it->second.position-genomic_position.position;

            {   // check the match between the SNV context and the corresponding nucleic 
                // triplet in `nucleic_sequence`

                // identify the first position in `nucleic_sequence` of the context
                size_t first_nucleotide = (SNV_pos==0?0:SNV_pos-1);

                // extract the `nucleic_sequence` nucleotides corresponding to the context
                std::string real_context = nucleic_sequence.substr(first_nucleotide, (SNV_pos==0?2:3));
                for (auto str_c: real_context) {
                    str_c = toupper(str_c);
                }

                // extract the part of the SNV context falling inside ``nucleic_sequence` 
                size_t length = real_context.size();
                std::string context = snv_it->second.context.get_sequence().substr((SNV_pos==0?1:0), length);

                // if they differ
                if (real_context!=context) {
                    std::ostringstream oss;

                    oss << "Context mismatch in chr. "<< GenomicPosition::chrtos(genomic_position.chr_id) 
                        << " in position " << snv_it->second.position << ": \"" 
                        << context << "\" expected and got \"" 
                        << real_context << "\"."<< std::endl 
                        << "Are you using a FASTA file different from that used for the context index?"
                        << std::endl;
                    throw std::domain_error(oss.str());
                }
            }

            // change the `nucleic_sequence` nucleotide according to the considered SNV
            nucleic_sequence[SNV_pos] = snv_it->second.mutated_base;
        }

        return nucleic_sequence;
    }

    /**
     * @brief Read the next chromosome data from the FASTA stream
     * 
     * @tparam DATA_TYPE if the type of the data to be read from the stream
     * @param FASTA_stream is the FASTA stream
     * @param chr_data is the object that will be filled with the chromosome data if
     *          some chromosome is read from `FASTA_stream`
     * @param progress_bar is the progress bar 
     * @return `true` if and only if a chromosome sequence is read from `FASTA_stream`
     */
    template<typename DATA_TYPE, std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
    static bool read_next_chromosome_data(std::istream& FASTA_stream, ChromosomeData<DATA_TYPE>& chr_data, 
                                          UI::ProgressBar& progress_bar)
    {
        using namespace Races::IO::FASTA;

        FilterNonChromosomeSequence filter;

        if (DATA_TYPE::read(FASTA_stream, chr_data, filter, progress_bar)) {
            chr_data.chr_id = filter.last_chr_id;

            return true;
        }

        return false;
    }

    /**
     * @brief Write the SAM alignment line corresponding to a simulated read
     * 
     * This method writes in a stream the SAM alignment line corresponding to a simulated read
     * whose initial position on a reference genome chromosome, and SNVs are provided.
     * 
     * @tparam DATA_TYPE is the type of information available about the considered chromosome
     * @param SAM_stream is the SAM file output stream
     * @param chr_data is the data about the chromosome from which the simulated read come from
     * @param SNVs is a map from genomic positions to the corresponding SNVs
     * @param read_first_position is the first position of the simulated read
     */
    template<typename DATA_TYPE, 
         std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
    void write_SAM_alignment(std::ostream& SAM_stream, const ChromosomeData<DATA_TYPE>& chr_data,
                             const std::map<GenomicPosition, SNV>& SNVs,
                             const ChrPosition& read_first_position)
    {
        const GenomicPosition genomic_position{chr_data.chr_id, read_first_position};

        std::string read;
        if constexpr(std::is_base_of_v<Races::IO::FASTA::Sequence, DATA_TYPE>) {
            read = chr_data.nucleotides.substr(read_first_position-1, read_size);
            apply_SNVs(read, SNVs, genomic_position);
        } else {
            read = "*";
        }

        const auto mapq = 33;

        SAM_stream << "r" <<  std::setfill('0') << std::setw(11) << ++num_of_reads; // QNAME
        SAM_stream << '\t' << '0'                  // FLAG
                   << '\t' << chr_data.name       // RNAME
                   << '\t' << read_first_position  // POS
                   << '\t' << mapq                 // MAPQ
                   << '\t' << get_CIGAR(SNVs, genomic_position, read_size) // CIGAR
                   << '\t' << '*'   // RNEXT
                   << '\t' << '0'   // PNEXT
                   << '\t' << '0'   // TLEN
                   << '\t' << read  // SEQ
                   << '\t' << '*'   // QUAL
                   << std::endl;
    }

    /**
     * @brief Generate simulated reads on a chromosome and write their SAM alignments
     * 
     * This method takes a list of genome mutations and a reference genome chromosome, it 
     * generates simulated reads over it up to a specified coverage of the mutated genomes, 
     * and write the corresponding SAM alignments in a stream.
     * 
     * @tparam DATA_TYPE is the type of information available about the considered chromosome
     * @param SAM_stream is the SAM file output stream
     * @param mutations is a list of genome mutations
     * @param chr_data is the data about the chromosome from which the simulated read come from
     * @param allelic_coverage is the aimed coverage the alleles of the mutated genomes
     * @param total_steps is the total number of steps required to complete the overall procedure
     * @param steps is the number of performed steps
     * @param progress_bar is the progress bar 
     */
    template<typename DATA_TYPE, std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
    void generate_SAM_alignments(std::ostream& SAM_stream, const std::list<GenomeMutations>& mutations, 
                                 const ChromosomeData<DATA_TYPE>& chr_data,
                                 const double& allelic_coverage, const size_t& total_steps, 
                                 size_t& steps, Races::UI::ProgressBar& progress_bar)
    {
        progress_bar.set_message("Processing chr. "+GenomicPosition::chrtos(chr_data.chr_id));

        for (const auto& cell_mutations: mutations) {
            const auto& chr_mutations = cell_mutations.get_chromosome(chr_data.chr_id);

            for (const auto& [allele_id, allele] : chr_mutations.get_alleles()) {
                for (const auto& [position, fragment] : allele.get_fragments()) {
                    if (fragment.get_final_position()>read_size) {
                        std::poisson_distribution<size_t> p_dist(allelic_coverage*fragment.size()/read_size);

                        size_t num_of_reads = p_dist(random_generator);

                        const auto& SNVs = fragment.get_SNVs();

                        auto first_possible_begin = fragment.get_initial_position();
                        auto last_possible_begin = static_cast<ChrPosition>(fragment.get_final_position()-read_size+1);

                        std::uniform_int_distribution<ChrPosition> dist(first_possible_begin, last_possible_begin);
                        for (size_t i=0; i<num_of_reads; ++i) {
                            write_SAM_alignment(SAM_stream, chr_data, SNVs, dist(random_generator));
                        }
                    }
                }

                steps += allele.size();

                progress_bar.set_progress(100*steps/total_steps);
            }
        }
    }

    /**
     * @brief Generate simulated reads and write their SAM alignments
     * 
     * This method takes a list of genome mutations, it generates simulated reads up to a 
     * specified coverage of the mutated genomes, and write the corresponding SAM alignments
     * in files, named after genome chromosomes, inside an output directory. 
     * 
     * @tparam DATA_TYPE is the type of information available about the considered chromosome
     * @param SAM_directory_name is the name of the directory for the output SAM files 
     * @param mutations is a list of genome mutations
     * @param allelic_coverage is the aimed coverage the alleles of the mutated genomes
     * @param total_steps is the total number of steps required to complete the overall procedure
     * @param steps is the number of performed steps
     * @param progress_bar is the progress bar
     */
    template<typename DATA_TYPE=Races::IO::FASTA::Sequence,
             std::enable_if_t<std::is_base_of_v<Races::IO::FASTA::SequenceInfo, DATA_TYPE>, bool> = true>
    void generate_SAM(const std::list<GenomeMutations>& mutations, const double& allelic_coverage,
                      const size_t& total_steps, size_t& steps, UI::ProgressBar& progress_bar)
    {
        std::ifstream ref_stream(ref_genome_filename);

        ChromosomeData<DATA_TYPE> chr_data;
        
        progress_bar.set_progress((100*steps)/total_steps, "Reading next chromosome");
        while (read_next_chromosome_data(ref_stream, chr_data, progress_bar)) {
            progress_bar.set_progress((100*(++steps))/total_steps);

            auto chr_str = GenomicPosition::chrtos(chr_data.chr_id);

            std::string SAM_filename = "chr_"+chr_str+".sam";
            auto SAM_file_path = SAM_directory_name/SAM_filename;
            std::ofstream SAM_stream;

            if (!std::filesystem::exists(SAM_file_path)) {
                SAM_stream = std::ofstream(SAM_file_path);

                SAM_stream << "@HD\tVN:1.6\tSO:unknown" << std::endl;
                SAM_stream << "@SQ\tSN:" << chr_data.name << "\tLN:" << chr_data.length 
                        << "\tAN:chromosome" << chr_str << ",chr" << chr_str 
                        << ",chromosome_" << chr_str << ",chr_" << chr_str << std::endl;
            } else {
                if (!std::filesystem::is_regular_file(SAM_file_path)) {
                    throw std::runtime_error("\""+std::string(SAM_file_path)+"\" is not a regular file");
                }

                SAM_stream = std::ofstream(SAM_file_path, std::ofstream::app);
            }
            generate_SAM_alignments(SAM_stream, mutations, chr_data, allelic_coverage, 
                                    total_steps, steps, progress_bar);
        }
        
        progress_bar.set_message("SAM Produced");
        progress_bar.set_progress(100);
    }
public:
    /**
     * @brief SAM generator mode
     */
    enum class Mode {
        CREATE,     // Create
        OVERWRITE,  // Overwrite files
        APPEND      // Append files
    };

    /**
     * @brief The empty constructor
     */
    SAMGenerator():
        random_generator(), SAM_directory_name(""), ref_genome_filename(""), read_size(0), num_of_reads(0)
    {}

    /**
     * @brief Create a SAM generator
     * 
     * @param SAM_directory_name is the SAM output directory
     * @param ref_genome_filename is reference genome filename
     * @param read_size is the size of the output reads
     * @param seed is the random generator seed
     * @param mode is the SAM generator output mode
     */
    SAMGenerator(const std::string& SAM_directory_name, const std::string& ref_genome_filename,
                 const size_t& read_size, const int& seed=0, const Mode mode=Mode::CREATE):
        random_generator(seed), SAM_directory_name(SAM_directory_name), ref_genome_filename(ref_genome_filename), 
        read_size(read_size),  num_of_reads(0)
    {
        namespace fs = std::filesystem;

        if (fs::exists(this->SAM_directory_name)) {
            if (mode==Mode::CREATE) {
                throw std::domain_error("\""+SAM_directory_name+"\" already exists");
            }

            if (!fs::is_directory(this->SAM_directory_name)) {
                throw std::domain_error("\""+SAM_directory_name+"\" is not a directory");
            }

            if (mode==Mode::OVERWRITE) {
                fs::remove_all(this->SAM_directory_name);
            }
        }

        fs::create_directory(this->SAM_directory_name);

        if (!fs::exists(this->ref_genome_filename)) {
            throw std::domain_error("\""+ref_genome_filename+"\" does not exist");
        }

        if (!fs::is_regular_file(this->ref_genome_filename)) {
            throw std::domain_error("\""+ref_genome_filename+"\" is not a regular file");
        }

        std::ifstream ref_stream(ref_genome_filename);
        char c = ref_stream.get();

        if (c != EOF && c != '>') {
            throw std::domain_error("\""+ref_genome_filename+"\" is not a FASTA file");
        }
    }

    /**
     * @brief Generate simulated reads and write the corresponding SAM alignments
     * 
     * This method generates simulated reads for a set of cell mutated genomes and 
     * writes the corresponding SAM alignments. The number of simulated reads 
     * depends on the specified coverage. Upon request, the SAM alignments avoid 
     * to report the reads sequence.
     * 
     * @param mutations is a list of genome mutations
     * @param coverage is the aimed coverage
     * @param with_sequences is a Boolean flag to include/avoid read sequence in SAM files
     * @param quiet is a Boolean flag to avoid progress bar and user messages
     * @return a reference to the updated object
     */
    SAMGenerator& operator()(const std::list<GenomeMutations>& mutations, 
                             const double& coverage, const bool with_sequences=true,
                             const bool quiet=false)
    {
        using namespace Races;

        if (coverage < 0) {
            std::ostringstream oss;

            oss << "coverage must be non-negative. Got " << coverage;
            throw std::domain_error(oss.str());
        }

        if (mutations.size()==0) {
            return *this;
        }

        UI::ProgressBar progress_bar(quiet);

        if (!represents_same_genome(mutations, progress_bar)) {
            throw std::domain_error("Not all the mutations represent the same genome");
        }

        auto overall_allelic_size = get_overall_allelic_size(mutations, progress_bar);
        auto allelic_coverage = (coverage*mutations.front().size())/overall_allelic_size;

        auto total_steps = overall_allelic_size + 2*(mutations.front().get_chromosomes().size())*mutations.size();
        size_t steps{0};

        if (with_sequences) {
            generate_SAM<Races::IO::FASTA::Sequence>(mutations, allelic_coverage,
                                                     total_steps, steps, progress_bar);
        } else {
            generate_SAM<Races::IO::FASTA::SequenceInfo>(mutations, allelic_coverage,
                                                         total_steps, steps, progress_bar);
        }

        return *this;
    }
};

}   // IO

}

}   // Races

#endif // __RACES_SAM_GENERATOR__