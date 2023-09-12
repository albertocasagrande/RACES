/**
 * @file genome_mutations.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implements genome and chromosome data structures
 * @version 0.7
 * @date 2023-09-12
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
    identifier(0), length(0), allelic_length(0), alleles(), next_allele_id(0)
{}

ChromosomeMutations::ChromosomeMutations(const ChromosomeId& identifier, const Length& size, const size_t& num_of_alleles):
    identifier(identifier), length(size), allelic_length(size*num_of_alleles), alleles(), next_allele_id(0)
{
    for (next_allele_id=0; next_allele_id<num_of_alleles; ++next_allele_id) {
        alleles.insert({next_allele_id, Allele{identifier, 1, size}});
    }
}

ChromosomeMutations::ChromosomeMutations(const GenomicRegion& chromosome_region, const size_t& num_of_alleles):
    identifier(chromosome_region.get_chromosome_id()), length(chromosome_region.size()), 
    allelic_length(chromosome_region.size()*num_of_alleles), alleles(), next_allele_id(0)        
{
    for (next_allele_id=0; next_allele_id<num_of_alleles; ++next_allele_id) {
        alleles.insert({next_allele_id, Allele{chromosome_region}});
    }
}

std::set<AlleleId> ChromosomeMutations::get_alleles_containing(const GenomicPosition& genomic_position) const
{
    std::set<AlleleId> allele_ids;

    for (const auto& [allele_id, allele]: alleles) {
        if (allele.contains(genomic_position)) {
            allele_ids.insert(allele_id);
        }
    }

    return allele_ids;
}

std::set<AlleleId> ChromosomeMutations::get_alleles_containing(const GenomicRegion& genomic_region) const
{
    std::set<AlleleId> allele_ids;

    for (const auto& [allele_id, allele]: alleles) {
        if (allele.contains(genomic_region)) {
            allele_ids.insert(allele_id);
        }
    }

    return allele_ids;
}

const Allele& ChromosomeMutations::find_allele(const AlleleId& allele_id) const
{
    auto it = alleles.find(allele_id);

    if (it == alleles.end()) {
        std::ostringstream oss;

        oss << "Chromosome " << identifier <<  " has not allele " << allele_id << ".";
        throw std::out_of_range(oss.str());
    }

    return it->second;
}

Allele& ChromosomeMutations::find_allele(const AlleleId& allele_id)
{
    auto it = alleles.find(allele_id);

    if (it == alleles.end()) {
        std::ostringstream oss;

        oss << "Chromosome " << identifier <<  " has not allele " << allele_id << ".";
        throw std::out_of_range(oss.str());
    }

    return it->second;
}

bool ChromosomeMutations::allele_contains(const AlleleId& allele_id, const GenomicRegion& genomic_region) const
{
    if (!contains(genomic_region)) {
        throw std::domain_error("The genomic region is not in the chromosome");
    }

    return find_allele(allele_id).contains(genomic_region);
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

bool ChromosomeMutations::amplify_region(const GenomicRegion& genomic_region, const AlleleId& allele_id)
{
    if (!allele_contains(allele_id, genomic_region)) {
        return false;
    }

    auto new_allele = alleles.at(allele_id).copy(genomic_region);

    allelic_length += new_allele.size();

    alleles.insert({next_allele_id, std::move(new_allele)});

    CNAs.push_back(CopyNumberAlteration::new_amplification(genomic_region, allele_id, next_allele_id));

    ++next_allele_id;

    return true;
}


bool ChromosomeMutations::remove_region(const GenomicRegion& genomic_region, const AlleleId& allele_id)
{
    if (!allele_contains(allele_id, genomic_region)) {
        return false;
    }

    alleles[allele_id].remove(genomic_region);

    CNAs.push_back(CopyNumberAlteration::new_deletion(genomic_region, allele_id));

    return true;
}

bool ChromosomeMutations::has_context_free(const GenomicPosition& genomic_position) const
{
    if (!contains(genomic_position)) {
        throw std::domain_error("The genomic position is not in the chromosome");
    }

    for (const auto& [allele_id, allele]: alleles) {
        if (!allele.has_context_free(genomic_position)) {
            return false;
        }
    }

    return true;
}

bool ChromosomeMutations::insert(const SNV& snv, const AlleleId& allele_id)
{
    if (!contains(snv)) {
        throw std::domain_error("The genomic position of the SNV is not in the chromosome");
    }

    if (!has_context_free(snv)) {
        return false;
    }

    return find_allele(allele_id).insert(snv);
}

bool ChromosomeMutations::remove_SNV(const GenomicPosition& genomic_position)
{
    if (!contains(genomic_position)) {
        throw std::domain_error("The genomic position of the SNV is not in the chromosome");
    }

    for (auto& [allele_id, allele]: alleles) {
        if (allele.remove_SNV(genomic_position)) {
            return true;
        }
    }

    return false;
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

std::set<AlleleId> GenomeMutations::get_alleles_containing(const GenomicPosition& genomic_position) const
{
    auto chr_it = find_chromosome(chromosomes, genomic_position.chr_id);

    return chr_it->second.get_alleles_containing(genomic_position);
}

std::set<AlleleId> GenomeMutations::get_alleles_containing(const GenomicRegion& genomic_region) const
{
    auto chr_it = find_chromosome(chromosomes, genomic_region.get_chromosome_id());

    return chr_it->second.get_alleles_containing(genomic_region);
}

bool GenomeMutations::allele_contains(const AlleleId& allele_id, const GenomicRegion& genomic_region) const
{
    auto chr_it = find_chromosome(chromosomes, genomic_region.get_chromosome_id());

    return chr_it->second.allele_contains(allele_id, genomic_region);
}

bool GenomeMutations::amplify_region(const GenomicRegion& genomic_region, const AlleleId& allele_id)
{
    auto chr_it = find_chromosome(chromosomes, genomic_region.get_chromosome_id());

    return chr_it->second.amplify_region(genomic_region, allele_id);
}

bool GenomeMutations::remove_region(const GenomicRegion& genomic_region, const AlleleId& allele_id)
{
    auto chr_it = find_chromosome(chromosomes, genomic_region.get_chromosome_id());

    return chr_it->second.remove_region(genomic_region, allele_id);
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

CellGenomeMutations::CellGenomeMutations():
    Cell(), GenomeMutations()
{}


CellGenomeMutations::CellGenomeMutations(const Drivers::Cell& cell, const GenomeMutations& genome_mutations):
    Cell(cell), GenomeMutations(genome_mutations)
{}

}   // Passengers

}   // Races

namespace std
{

std::ostream& operator<<(std::ostream& os, const Races::Passengers::ChromosomeMutations& chromosome_mutations)
{
    os << "Chromosome " << Races::Passengers::GenomicPosition::chrtos(chromosome_mutations.id()) << std::endl;

    for (const auto& [allele_id, allele]: chromosome_mutations.get_alleles()) {
        os << "  " << allele_id << ": " << allele << std::endl;
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const Races::Passengers::GenomeMutations& genome_mutations)
{
    for (const auto& [chr_id, chr_mutations]: genome_mutations.get_chromosomes()) {
        os << chr_mutations << std::endl;
    }

    return os;
}

}   // std