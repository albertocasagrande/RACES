/**
 * @file genome_mutations.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implements genome and chromosome data structures
 * @version 0.1
 * @date 2023-07-27
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

#include <sstream>
#include <cassert>

#include "genome_mutations.hpp"

namespace Races 
{

namespace Passengers
{

ChromosomeMutations::ChromosomeMutations():
    identifier(0), length(0), allelic_length(0), fragments(), next_allele_id(0)
{}

ChromosomeMutations::ChromosomeMutations(const ChromosomeId& identifier, const Length& size, const size_t& num_of_alleles):
    identifier(identifier), length(size), allelic_length(size*num_of_alleles), fragments(), next_allele_id(num_of_alleles)
{
    Fragment chr_fragment{identifier, size, num_of_alleles};

    fragments[chr_fragment.get_begin()] = chr_fragment;
}

ChromosomeMutations::ChromosomeMutations(const GenomicRegion& chromosome_region, const size_t& num_of_alleles):
    identifier(chromosome_region.get_chromosome_id()), length(chromosome_region.size()), 
    allelic_length(chromosome_region.size()*num_of_alleles), fragments(), next_allele_id(num_of_alleles)        
{
    fragments[chromosome_region.get_begin()] = Fragment(chromosome_region, num_of_alleles);
}

auto find_last_not_starting_after(std::map<GenomicPosition, Fragment>& fragments, const GenomicPosition& genomic_position)
{
    // find the first fragment that begins *after* `genomic_position`
    auto it = fragments.upper_bound(genomic_position);

    assert(it == fragments.end() || it->first.position > genomic_position.position);

    // if its not the first in the chromosome
    if (it != fragments.begin()) {

        // consider the previous one, i.e., the last one that started 
        // before `genomic_position`
        --it;
    }

    return it;
}

auto find_first_not_ending_before(std::map<GenomicPosition, Fragment>& fragments, const GenomicPosition& genomic_position)
{
    auto it = find_last_not_starting_after(fragments, genomic_position);

    // if it ends before `genomic_region`
    if (it != fragments.end() && it->second.ends_before(genomic_position)) {

        // consider the following one
        ++it;
    }

    return it;
}


auto find_last_not_starting_after(const std::map<GenomicPosition, Fragment>& fragments, const GenomicPosition& genomic_position)
{
    // find the first fragment that begins *after* `genomic_position`
    auto it = fragments.upper_bound(genomic_position);

    assert(it == fragments.end() || it->first.position > genomic_position.position);

    // if its not the first in the chromosome
    if (it != fragments.begin()) {

        // consider the previous one, i.e., the last one that started 
        // before `genomic_position`
        --it;
    }

    return it;
}

auto find_first_not_ending_before(const std::map<GenomicPosition, Fragment>& fragments, const GenomicPosition& genomic_position)
{
    auto it = find_last_not_starting_after(fragments, genomic_position);

    // if it ends before `genomic_region`
    if (it != fragments.end() && it->second.ends_before(genomic_position)) {

        // consider the following one
        ++it;
    }

    return it;
}

std::set<AlleleId> ChromosomeMutations::get_allele_ids_in(const GenomicRegion& genomic_region) const
{
    // find the first fragment that may overlaps genomic region
    auto f_it = find_first_not_ending_before(fragments, genomic_region.get_begin());
    
    if (f_it == fragments.end() || f_it->second.begins_after(genomic_region.get_begin())) {

        // no fragment touches the genomic_region
        return {};
    }

    // get the allele ids from the first overlapping 
    // fragment
    auto allele_ids = f_it->second.get_allele_ids();

    // for all the overlapping fragments
    for (; f_it!=fragments.end() && f_it->second.overlaps(genomic_region) 
            && !allele_ids.empty(); ++f_it) {

        // remove all the ids that are not contained 
        // by the fragment
        auto ids_it = allele_ids.begin();
        while (ids_it != allele_ids.end()) {
            if (!f_it->second.has_allele(*ids_it)) {
                allele_ids.extract(ids_it++);
            } else {
                ++ids_it;
            }
        }
    }

    return allele_ids;
}

std::map<AlleleId, std::set<SNV>> ChromosomeMutations::get_SNVs_per_allele() const
{
    std::map<AlleleId, std::set<SNV>> SNVs;

    for (const auto& [begin_pos, fragment]: fragments) {
        for (const auto& [allele_id, allele]: fragment.get_alleles()) {
            auto& allele_SNVs = SNVs[allele_id];
            for (const auto& [pos, snv]: allele) {
                allele_SNVs.insert(snv);
            }
        }
    }

    return SNVs;
}

void split_fragments_overlapping(std::map<GenomicPosition, Fragment>& fragments, const GenomicRegion& genomic_region)
{
    auto f_it = find_first_not_ending_before(fragments, genomic_region.get_begin());

    // if some fragment not ending before `genomic_region` begin exists and 
    // it does begin before the initial position of `genomic_region`
    if (f_it != fragments.end()
            && f_it->second.get_initial_position()<genomic_region.get_initial_position()) {
        
        // split it 
        auto new_fragment = f_it->second.split(genomic_region.get_begin());

        // add the new fragment to the chromosome fragments
        fragments[new_fragment.get_begin()] = new_fragment;
    }

    auto region_end = genomic_region.get_end();
    f_it = find_first_not_ending_before(fragments, region_end);

    // if some fragment ending after `genomic_region` end exists and 
    // it does begin before or at the final position of `genomic_region`
    if (f_it != fragments.end() 
            && f_it->second.get_initial_position()<=genomic_region.get_final_position()
            && f_it->second.get_final_position()>genomic_region.get_final_position()) {
        
        // build the split point one nucleotide after the end of `genomic_region`
        GenomicPosition split_point(region_end.chr_id, region_end.position+1);

        // split the fragment
        auto new_fragment = f_it->second.split(split_point);

        // add the new fragment to the chromosome fragments
        fragments[new_fragment.get_begin()] = new_fragment;
    }
}

bool ChromosomeMutations::contains(const GenomicPosition& genomic_position) const
{
    return (genomic_position.chr_id==id()       // same chromosome
            && genomic_position.position>0      // lays in [1, chromosome size]
            && genomic_position.position<=size());
}

bool ChromosomeMutations::contains(const GenomicRegion& genomic_region) const
{
    return (genomic_region.get_chromosome_id()==id()     // same chromosome
            && genomic_region.get_initial_position()>0   // the lays in [1, chromosome size]
            && genomic_region.get_final_position()<=size());
}

bool ChromosomeMutations::has_allele_on(const AlleleId& allele_id, const GenomicRegion& genomic_region) const
{
    if (!contains(genomic_region)) {
        throw std::domain_error("The genomic region is not in the chromosome");
    }

    auto f_it = find_first_not_ending_before(fragments, genomic_region.get_begin());

    assert(f_it->second.get_final_position() >= genomic_region.get_initial_position());

    if (f_it == fragments.end() || 
        f_it->second.get_initial_position() > genomic_region.get_initial_position()) {
        return false;
    }

    // check that all the overlapping fragments contains allele_ids and 
    // that they are continuous
    auto prec_it=f_it;
    while (f_it!=fragments.end() && f_it->second.overlaps(genomic_region)) {
        if (!f_it->second.has_allele(allele_id)
            ||(prec_it!=f_it 
               && (prec_it->second.get_final_position()+1 < f_it->second.get_initial_position()))) {
            return false;
        }

        prec_it=f_it;
        ++f_it;
    }

    return true;
}

bool ChromosomeMutations::increase_copy_number(const GenomicRegion& genomic_region, const AlleleId& allele_id)
{
    if (!has_allele_on(allele_id, genomic_region)) {
        return false;
    }

    split_fragments_overlapping(fragments, genomic_region);

    // duplicate the allele over all the fragments
    for (auto f_it = find_first_not_ending_before(fragments, genomic_region.get_begin());
            f_it!=fragments.end() && !genomic_region.ends_before(f_it->first); ++f_it) {

        f_it->second.duplicate_allele(allele_id, next_allele_id);
    }
    ++next_allele_id;

    allelic_length += genomic_region.size();

    return true;
}


bool ChromosomeMutations::decrease_copy_number(const GenomicRegion& genomic_region, const AlleleId& allele_id)
{
    if (!has_allele_on(allele_id, genomic_region)) {
        return false;
    }

    // remove the specified allele
    for (auto f_it = find_first_not_ending_before(fragments, genomic_region.get_begin());
            f_it!=fragments.end() && !genomic_region.ends_before(f_it->first); ++f_it) {
        
        f_it->second.remove_allele(allele_id);
    }
    
    allelic_length -= genomic_region.size();

    return true;
}

bool ChromosomeMutations::is_mutational_context_free(GenomicPosition genomic_position) const
{
    if (genomic_position.position < 2) {
        return false;
    }
 
    genomic_position.position -= 1;

    auto f_it = find_first_not_ending_before(fragments, genomic_position);

    for (size_t i=0; i<3; ++i) {
        genomic_position.position += i;
        if (f_it->second.get_final_position()<genomic_position.position) {
            ++f_it;

            if (f_it==fragments.end()) {
                return false;
            }
        }

        if ((f_it!=fragments.end() && f_it->second.get_initial_position()>genomic_position.position)
                || (f_it->second.has_SNV_at(genomic_position))) {
            return false;
        }
    }

    return true;
}

bool ChromosomeMutations::insert(const SNV& snv, const AlleleId allele_id)
{
    if (!contains(snv)) {
        throw std::domain_error("The genomic position of the SNV is not in the chromosome");
    }

    if (!is_mutational_context_free(snv)) {
        std::cout << "QUI " << snv << std::endl;
        return false;
    }

    auto f_it = find_first_not_ending_before(fragments, snv);
    if (f_it->first.position <= snv.position) {
        return f_it->second.insert(snv, allele_id);   
    }

    return false;
}

bool ChromosomeMutations::remove_SNV(const GenomicPosition& genomic_position)
{
    if (!contains(genomic_position)) {
        throw std::domain_error("The genomic position of the SNV is not in the chromosome");
    }

    auto f_it = find_first_not_ending_before(fragments, genomic_position);

    if (f_it->first.position > genomic_position.position) {
        return false;
    }

    return f_it->second.remove_SNV(genomic_position);
}

GenomeMutations::GenomeMutations()
{}

GenomeMutations::GenomeMutations(const std::vector<ChromosomeMutations>& chromosomes)
{
    for (const auto& chromosome: chromosomes) {
        if (this->chromosomes.count(chromosome.id())>0) {
            std::ostringstream oss;

            oss << "Chromosome " << GenomicPosition::chrtos(chromosome.id())
                << " has been specified twice"; 

            throw std::domain_error(oss.str());
        }
        this->chromosomes[chromosome.id()] = chromosome;
    }
}

GenomeMutations::Length GenomeMutations::size() const
{
    Length size = 0;

    for (const auto& [chr_id, chromosome]: chromosomes) {
        size += chromosome.size();
    }

    return size;
}

GenomeMutations::Length GenomeMutations::allelic_size() const
{
    Length allelic_size = 0;

    for (const auto& [chr_id, chromosome]: chromosomes) {
        allelic_size += chromosome.allelic_size();
    }

    return allelic_size;
}

auto find_chromosome(std::map<ChromosomeId, ChromosomeMutations>& chromosomes, const ChromosomeId& chr_id)
{
    auto chr_it = chromosomes.find(chr_id);

    if (chr_it == chromosomes.end()) {
        std::ostringstream oss;

        oss << "Unknown chromosome " << GenomicPosition::chrtos(chr_id); 
        throw std::domain_error(oss.str());
    }

    return chr_it;   
}

auto find_chromosome(const std::map<ChromosomeId, ChromosomeMutations>& chromosomes, const ChromosomeId& chr_id)
{
    auto chr_it = chromosomes.find(chr_id);

    if (chr_it == chromosomes.end()) {
        std::ostringstream oss;

        oss << "Unknown chromosome " << GenomicPosition::chrtos(chr_id); 
        throw std::domain_error(oss.str());
    }

    return chr_it;   
}

std::set<AlleleId> GenomeMutations::get_allele_ids_in(const GenomicRegion& genomic_region) const
{
    auto chr_it = find_chromosome(chromosomes, genomic_region.get_chromosome_id());

    return chr_it->second.get_allele_ids_in(genomic_region);
}

bool GenomeMutations::has_allele_on(const AlleleId& allele_id, const GenomicRegion& genomic_region) const
{
    auto chr_it = find_chromosome(chromosomes, genomic_region.get_chromosome_id());

    return chr_it->second.has_allele_on(allele_id, genomic_region);
}

bool GenomeMutations::increase_copy_number(const GenomicRegion& genomic_region, const AlleleId& allele_id)
{
    auto chr_it = find_chromosome(chromosomes, genomic_region.get_chromosome_id());

    return chr_it->second.increase_copy_number(genomic_region, allele_id);
}

bool GenomeMutations::decrease_copy_number(const GenomicRegion& genomic_region, const AlleleId& allele_id)
{
    auto chr_it = find_chromosome(chromosomes, genomic_region.get_chromosome_id());

    return chr_it->second.decrease_copy_number(genomic_region, allele_id);
}

bool GenomeMutations::insert(const SNV& snv, const AlleleId& allele_id)
{
    auto chr_it = find_chromosome(chromosomes, snv.chr_id);

    return chr_it->second.insert(snv, allele_id);
}

bool GenomeMutations::remove_SNV(const GenomicPosition& genomic_position)
{
    auto chr_it = find_chromosome(chromosomes, genomic_position.chr_id);

    return chr_it->second.remove_SNV(genomic_position);
}

}   // Passengers

}   // Races

