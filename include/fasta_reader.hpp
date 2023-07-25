/**
 * @file fasta_reader.hpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Defines a FASTA file reader and support structures
 * @version 0.2
 * @date 2023-07-25
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

#ifndef __RACES_FASTA_READER__
#define __RACES_FASTA_READER__

#include <cstdint>
#include <string>
#include <istream>

#include "fragment.hpp"

namespace Races
{

namespace IO
{

namespace FASTA
{


/**
 * @brief A structure to decode human genome chromosome sequences name 
 */
struct HumanSeqNameDecoder
{
    /**
     * @brief Turn the human chromosome name from string for to chromosome id
     * 
     * This method performs a standard transformation from string to integer
     * for all the human chromosomes, but 'X' and 'Y' whose chromosome ids 
     * are 23 and 24, respectively.
     * 
     * @param chr_name is a string containing the human chromosome name
     * @return the chromosome id associated to `chr_name`
     */
    static Passengers::ChromosomeId stochr(const std::string& chr_name);

    /**
     * @brief Turn a human chromosome id into the corresponding chromosome name
     * 
     * This method performs a standard transformation from integer to string
     * for all the human chromosomes, but 'X' and 'Y' whose chromosome ids 
     * are 23 and 24, respectively.
     * 
     * @param chr_id is the chromosome id
     * @return the chromosome nome associated to `chr_id` 
     */
    static std::string chrtos(const Passengers::ChromosomeId& chr_id);

    /**
     * @brief Extract the genomic region from a sequence name
     * 
     * @param seq_name is a FASTA sequence name
     * @param chr_region is the variable where the chromosome region will be placed
     * @return `true` if and only if `seq_name` correspond to a DNA chromosome
     */
    static bool get_chromosome_region(const std::string& seq_name, Passengers::GenomicRegion& chr_region);
};

/**
 * @brief A class to read a sequence in a fasta file
 */
class SequenceReader
{
    std::istream *in;       //!< the reader input stream

    std::string seq_name;   //!< the sequence name
    uint32_t seq_pos;       //!< the sequence position

    unsigned int buffer_pos;    //!< the buffer position
    std::string buffer;         //!< the reader buffer
    bool concluded;             //!< the end-of-sequence flag

    /**
     * @brief Refill the buffer if necessary
     * 
     * This methods checks  whether the sequence has been fully read 
     * and, if this is not the case, it refills the buffer when 
     * necessary.
     * 
     * @return `true` if and only if the buffer has required refill
     */
    bool refill_buffer();
public:
    /**
     * @brief A sequence nucleotide input iterator
     */
    class nucleotide_iterator
    {
        SequenceReader* reader;  //!< the reader
        
        const bool valid_info;  //!< a flag to validated the sequence information
        std::string seq_name;   //!< the sequence name
        uint32_t seq_pos;       //!< the sequence position

        /**
         * @brief A constructor
         * 
         * @param reader is the sequence reader that build this
         *          nucleotide iterator
         */
        nucleotide_iterator(SequenceReader& reader);

    public:
        using difference_type   =   std::ptrdiff_t;
        using value_type        =   char;
        using pointer           =   const char*;
        using reference         =   const char&;
        using istream_type      =   std::istream;

        /**
         * @brief The empty constructor
         */
        nucleotide_iterator();

        /**
         * @brief Reference operator
         * 
         * @return a reference to the species pointer by the iterator 
         */
        reference operator*() const;

        /**
         * @brief Pointer operator
         * 
         * @return a pointer to the species pointer by the iterator 
         */
        pointer operator->() const;

        /**
         * @brief The prefix increment
         * 
         * @return a reference to the updated object
         */
        nucleotide_iterator& operator++();

        /**
         * @brief Get the name of the sequence been read
         * 
         * @return the name of the sequence
         */
        const std::string& get_sequence_name() const;

        /**
         * @brief Get the position of the read nucleotide 
         * 
         * @return the position of the read nucleotide  
         */
        const uint32_t& get_position() const;

        /**
         * @brief Test whether two iterator are pointing to the same nucleotide
         * 
         * @param it is nucleotide iterator
         * @return `true` if and only if both the iterator have reached the end 
         *      of the sequence of their are pointing to the same position in 
         *      the stream
         */
        inline bool operator==(const nucleotide_iterator& it) const
        {
            return ((reader==nullptr&&it.reader==nullptr)
                    || ((reader!=nullptr&&it.reader!=nullptr) &&
                        (reader->in==it.reader->in) && (seq_pos==it.seq_pos)));
        }

        /**
         * @brief Test whether two iterator are pointing to different nucleotide
         * 
         * @param it is nucleotide iterator
         * @return `false` if and only if both the iterator have reached the end 
         *      of the sequence of their are pointing to the same position in 
         *      the stream
         */
        inline bool operator!=(const nucleotide_iterator& it) const
        {
            return !(*this==it);
        }
    
        friend class SequenceReader;
    };

    /**
     * @brief The constructor
     * 
     * @param in is the input stream
     */
    SequenceReader(std::istream& in);

    /**
     * @brief Get the initial nucleotide iterator
     * 
     * @return the initial nucleotide iterator
     */
    inline nucleotide_iterator begin() {
        return nucleotide_iterator(*this);
    }

    /**
     * @brief Get the final nucleotide iterator
     * 
     * @return the final nucleotide iterator
     */
    inline nucleotide_iterator end() const {
        return nucleotide_iterator();
    }

    /**
     * @brief Remove the sequence from the input stream
     */
    void skip();

    /**
     * @brief Check whether the sequence has been fully read
     * 
     * This methods checks  whether the sequence has been fully read 
     * and, if this is not the case, it refills the buffer.
     * 
     * @return `true` if and only if the sequence has been fully read
     */
    inline bool finished() const
    {
        return concluded;
    }

    /**
     * @brief Get the position of the read nucleotide 
     * 
     * @return the position of the read nucleotide  
     */
    inline const std::string& get_sequence_name() const
    {
        return seq_name;
    }

    /**
     * @brief Get the name of the sequence been read
     * 
     * @return a constant reference to the name of the 
     *      sequence been read
     */
    inline const uint32_t& get_position() const
    {
        return seq_pos;
    }

    friend class nucleotide_iterator;
};

/**
 * @brief A class to read FASTA files
 */
class Reader
{
    std::istream& in;      //!< the input stream

public:
    /**
     * @brief A constructor
     * 
     * @param in is the input stream
     */
    Reader(std::istream& in);

    /**
     * @brief Test whether the FASTA file has been fully read
     * 
     * @return `true` no further sequences can be read from the 
     *          input stream
     */
    inline bool finished() {
        return in.eof();
    }

    /**
     * @brief Get a sequence reader for the next sequence
     * 
     * @return a sequence reader for the next sequence
     */
    SequenceReader next_sequence_reader();

    std::streampos tellg() const
    {
        return in.tellg();
    }
};

}   // FASTA

}   // IO

}   // Races

#endif // __RACES_FASTA_READER__