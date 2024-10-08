/**
 * @file genome_mutations.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements genome and chromosome data structures
 * @version 1.7
 * @date 2024-08-16
 *
 * @copyright Copyright (c) 2023-2024
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
#include <exception>

#include "genome_mutations.hpp"

namespace RACES
{

namespace Mutations
{

ChromosomeMutations::ChromosomeMutations():
    identifier(0), length(0), allelic_length(0), alleles(), next_allele_id(0)
{}

ChromosomeMutations::ChromosomeMutations(const ChromosomeId& identifier, const Length& size, const size_t& num_of_alleles):
    identifier(identifier), length(size), allelic_length(size*num_of_alleles), alleles(), next_allele_id(0)
{
    for (next_allele_id=0; next_allele_id<num_of_alleles; ++next_allele_id) {
        alleles.insert({next_allele_id, Allele{next_allele_id, identifier, 1, size}});
    }
}

ChromosomeMutations::ChromosomeMutations(const GenomicRegion& chromosome_region, const size_t& num_of_alleles):
    identifier(chromosome_region.get_chromosome_id()), length(chromosome_region.size()),
    allelic_length(chromosome_region.size()*num_of_alleles), alleles(), next_allele_id(0)
{
    for (next_allele_id=0; next_allele_id<num_of_alleles; ++next_allele_id) {
        alleles.insert({next_allele_id, Allele{next_allele_id, chromosome_region}});
    }
}

std::list<AlleleId> ChromosomeMutations::get_alleles_containing(const GenomicPosition& genomic_position) const
{
    std::list<AlleleId> allele_ids;

    for (const auto& [allele_id, allele]: alleles) {
        if (allele.contains(genomic_position)) {
            allele_ids.push_back(allele_id);
        }
    }

    return allele_ids;
}

std::list<AlleleId> ChromosomeMutations::get_alleles_containing(const GenomicRegion& genomic_region) const
{
    std::list<AlleleId> allele_ids;

    for (const auto& [allele_id, allele]: alleles) {
        if (allele.contains(genomic_region)) {
            allele_ids.push_back(allele_id);
        }
    }

    return allele_ids;
}

std::list<AlleleId> ChromosomeMutations::get_alleles_with_context_free_for(const SID& mutation) const
{
    std::list<AlleleId> context_free_alleles;

    for (const auto& [allele_id, allele]: alleles) {
        if (allele.has_context_free(mutation)) {
            context_free_alleles.push_back(allele_id);
        }
    }

    return context_free_alleles;
}

const Allele& ChromosomeMutations::get_allele(const AlleleId& allele_id) const
{
    auto it = alleles.find(allele_id);

    if (it == alleles.end()) {
        std::ostringstream oss;

        oss << "Chromosome " << GenomicPosition::chrtos(identifier)
            <<  " has not allele " << Allele::format_id(allele_id) << ".";
        throw std::out_of_range(oss.str());
    }

    return it->second;
}

Allele& ChromosomeMutations::get_allele(const AlleleId& allele_id)
{
    auto it = alleles.find(allele_id);

    if (it == alleles.end()) {
        std::ostringstream oss;

        oss << "Chromosome " << GenomicPosition::chrtos(identifier)
            <<  " has not allele " << Allele::format_id(allele_id) << ".";
        throw std::out_of_range(oss.str());
    }

    return it->second;
}

bool ChromosomeMutations::allele_contains(const AlleleId& allele_id, const GenomicRegion& genomic_region) const
{
    if (!contains(genomic_region)) {
        throw std::domain_error("The genomic region is not in the chromosome");
    }

    return get_allele(allele_id).contains(genomic_region);
}

bool ChromosomeMutations::contains(const SID& mutation) const
{
    const auto mutation_region = mutation.get_region();

    return contains(mutation_region);
}

bool ChromosomeMutations::amplify_region(const GenomicRegion& genomic_region, const AlleleId& allele_id,
                                         AlleleId& new_allele_id, const Mutation::Nature& nature)
{
    if (!contains(genomic_region)) {
        throw std::domain_error("The genomic region is not in the chromosome");
    }

    if (new_allele_id != RANDOM_ALLELE) {
        if (alleles.find(new_allele_id) != alleles.end()) {
            throw std::domain_error("Chromosome " + std::to_string(identifier)
                                    + " already contains the allele "
                                    + std::to_string(new_allele_id) + ".");
        }

        if (new_allele_id > next_allele_id) {
            next_allele_id = new_allele_id + 1;
        } else if (new_allele_id == next_allele_id) {
            ++(next_allele_id);
        }
    } else {
        new_allele_id = next_allele_id;

        ++(next_allele_id);
    }

    while (alleles.find(next_allele_id) != alleles.end()) {
        ++(next_allele_id);
    }

    const auto& allele = get_allele(allele_id);
    if (!allele.contains(genomic_region)) {
        return false;
    }

    auto new_allele = allele.copy(new_allele_id, genomic_region);

    allelic_length += new_allele.size();

    alleles.insert({new_allele_id, std::move(new_allele)});

    CNAs.push_back(CNA::new_amplification(genomic_region.get_begin(), genomic_region.size(),
                                          allele_id, new_allele_id, nature));

    return true;
}

bool ChromosomeMutations::remove_region(const GenomicRegion& genomic_region, const AlleleId& allele_id,
                                        const Mutation::Nature& nature)
{
    if (!contains(genomic_region)) {
        throw std::domain_error("The genomic region is not in the chromosome");
    }

    auto& allele = get_allele(allele_id);
    if (!allele.contains(genomic_region)) {
        return false;
    }

    allelic_length -= genomic_region.size();

    allele.remove(genomic_region);

    CNAs.push_back(CNA::new_deletion(genomic_region.get_begin(), genomic_region.size(),
                                     allele_id, nature));

    return true;
}

bool ChromosomeMutations::has_context_free(const SID& mutation) const
{
    if (!contains(mutation)) {
        std::ostringstream oss;

        oss << mutation << " is does not lays in the chromsome "
            << id() << ".";
        throw std::domain_error(oss.str());
    }

    bool some_allele_contains_pos{false};

    for (const auto& [allele_id, allele]: alleles) {
        if (allele.strictly_contains(mutation)) {
            some_allele_contains_pos = true;
            if (!allele.has_context_free(mutation)) {
                return false;
            }
        }
    }

    return some_allele_contains_pos;
}

bool ChromosomeMutations::apply_CNA(CNA& cna)
{
    const auto cna_region = cna.get_region();
    switch(cna.type) {
        case CNA::Type::AMPLIFICATION:
            return amplify_region(cna_region, cna.source,
                                  cna.dest, cna.nature);
        case CNA::Type::DELETION:
            return remove_region(cna_region, cna.dest,
                                 cna.nature);
        default:
            throw std::runtime_error("Unsupported CNA type");
    }
}

bool ChromosomeMutations::apply(const SID& mutation, const AlleleId& allele_id)
{
    if (!contains(mutation)) {
        throw std::domain_error("The genomic position of the SID is not in the chromosome");
    }

    Allele& allele = get_allele(allele_id);

    return allele.apply(mutation);
}

bool ChromosomeMutations::apply(const MutationSpec<SID>& mutation_spec)
{
    if (!contains(mutation_spec)) {
        throw std::domain_error("The genomic position of the SID is not in the chromosome");
    }

    Allele& allele = get_allele(mutation_spec.allele_id);

    return allele.apply(mutation_spec);
}

bool ChromosomeMutations::remove_mutation(const GenomicPosition& genomic_position)
{
    if (!contains(genomic_position)) {
        throw std::domain_error("The genomic position of the SID is not in the chromosome");
    }

    for (auto& [allele_id, allele]: alleles) {
        if (allele.remove_mutation(genomic_position)) {
            return true;
        }
    }

    return false;
}

bool ChromosomeMutations::includes(const SID& mutation) const
{
    if (mutation.chr_id != identifier) {
        return false;
    }

    for (auto& [allele_id, allele]: alleles) {
        if (allele.includes(mutation)) {
            return true;
        }
    }

    return false;
}

ChromosomeMutations ChromosomeMutations::copy_structure() const
{
    ChromosomeMutations copy;

    copy.identifier = identifier;
    copy.length = length;
    copy.allelic_length = allelic_length;
    copy.CNAs = CNAs;
    copy.next_allele_id = next_allele_id;

    for (const auto& [allele_id, allele] : alleles) {
        copy.alleles.emplace(allele_id, allele.copy_structure());
    }

    return copy;
}

void ChromosomeMutations::duplicate_alleles()
{
    auto it = alleles.begin();
    auto old_allele_id = next_allele_id;

    while (it != alleles.end()) {
        if (it->first < old_allele_id) {
            alleles.emplace(next_allele_id, it->second);
            ++next_allele_id;

            ++it;
        } else {
            it = alleles.end();
        }
    }

    allelic_length *= 2;
}

std::map<ChrPosition, AllelicType>
ChromosomeMutations::get_allelic_types(const std::set<ChrPosition>& break_points,
                                       const size_t& min_allelic_size) const
{
    size_t allelic_size{min_allelic_size};
    std::map<ChrPosition, std::map<AlleleId, uint32_t>> allele_counter;
    for (const auto& [allele_id, allele] : get_alleles()) {
        const auto orig_allele_id = allele.get_history().front();
        if (orig_allele_id+1>allelic_size) {
            allelic_size = orig_allele_id+1;
        }
        for (const auto& [pos, fragment] : allele.get_fragments()) {
            auto bp_it = break_points.find(pos.position);

            if (bp_it == break_points.end()) {
                std::ostringstream oss;

                oss << "The breakpoint " << pos << " is not in "
                    << "the provided breakpoint set.";
                throw std::domain_error(oss.str());
            }

            const auto fragment_end = fragment.get_final_position();
            while (*bp_it < fragment_end) {
                ++(allele_counter[*bp_it][orig_allele_id]);
                ++bp_it;
            }
        }
    }

    std::map<ChrPosition, AllelicType> allelic_map;
    for (const auto& [pos, counter] : allele_counter) {
        AllelicType atype(allelic_size, 0);

        for (const auto& [allele_id, count] : counter) {
            atype[allele_id] = count;
        }

        allelic_map[pos] = atype;
    }

    return allelic_map;
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

GenomeMutations::GenomeMutations(const std::list<GenomicRegion>& chromosome_regions,
                                 const size_t& num_of_alleles)
{
    for (const auto& chr_region : chromosome_regions) {
        chromosomes[chr_region.get_chromosome_id()] = ChromosomeMutations(chr_region, num_of_alleles);
    }
}

GenomeMutations::GenomeMutations(const std::list<GenomicRegion>& chromosome_regions,
                                 const std::map<ChromosomeId, size_t>& alleles_per_chromosome)
{
    for (const auto& chr_region : chromosome_regions) {
        auto chr_id = chr_region.get_chromosome_id();
        chromosomes[chr_id] = ChromosomeMutations(chr_region,
                                                    alleles_per_chromosome.at(chr_id));
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
    Length allelic_size{0};

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

std::map<ChromosomeId, size_t> GenomeMutations::get_absolute_chromosome_positions() const
{
    std::map<ChromosomeId, size_t> abs_pos;

    size_t current_pos = 1;
    for (const auto& [chr_id, chr]: chromosomes) {
        abs_pos[chr_id] = current_pos;
        current_pos += chr.size();
    }

    return abs_pos;
}

const ChromosomeMutations& GenomeMutations::get_chromosome(const ChromosomeId& chromosome_id) const
{
    auto found = chromosomes.find(chromosome_id);

    if (found == chromosomes.end()) {
        throw std::out_of_range("Unknown chromosome "
                                + GenomicPosition::chrtos(chromosome_id) + ".");
    }
    return found->second;
}

ChromosomeMutations& GenomeMutations::get_chromosome(const ChromosomeId& chromosome_id)
{
    auto found = chromosomes.find(chromosome_id);

    if (found == chromosomes.end()) {
        throw std::out_of_range("Unknown chromosome "
                                + GenomicPosition::chrtos(chromosome_id) + ".");
    }
    return found->second;
}

std::list<AlleleId> GenomeMutations::get_alleles_containing(const GenomicPosition& genomic_position) const
{
    auto chr_it = find_chromosome(chromosomes, genomic_position.chr_id);

    return chr_it->second.get_alleles_containing(genomic_position);
}

std::list<AlleleId> GenomeMutations::get_alleles_containing(const GenomicRegion& genomic_region) const
{
    auto chr_it = find_chromosome(chromosomes, genomic_region.get_chromosome_id());

    return chr_it->second.get_alleles_containing(genomic_region);
}

std::list<AlleleId> GenomeMutations::get_alleles_with_context_free_for(const SID& mutation) const
{
    auto chr_it = find_chromosome(chromosomes, mutation.chr_id);

    return chr_it->second.get_alleles_with_context_free_for(mutation);
}

bool GenomeMutations::allele_contains(const AlleleId& allele_id, const GenomicRegion& genomic_region) const
{
    auto chr_it = find_chromosome(chromosomes, genomic_region.get_chromosome_id());

    return chr_it->second.allele_contains(allele_id, genomic_region);
}

bool GenomeMutations::amplify_region(const GenomicRegion& genomic_region, const AlleleId& allele_id,
                                     AlleleId& new_allele_id, const Mutation::Nature& nature)
{
    auto chr_it = find_chromosome(chromosomes, genomic_region.get_chromosome_id());

    return chr_it->second.amplify_region(genomic_region, allele_id, new_allele_id, nature);
}

bool GenomeMutations::remove_region(const GenomicRegion& genomic_region, const AlleleId& allele_id,
                                    const Mutation::Nature& nature)
{
    auto chr_it = find_chromosome(chromosomes, genomic_region.get_chromosome_id());

    return chr_it->second.remove_region(genomic_region, allele_id, nature);
}

bool GenomeMutations::apply_CNA(CNA& cna)
{
    const auto cna_region = cna.get_region();
    switch(cna.type) {
        case CNA::Type::AMPLIFICATION:
            return amplify_region(cna_region, cna.source,
                                  cna.dest, cna.nature);
        case CNA::Type::DELETION:
            return remove_region(cna_region, cna.dest,
                                 cna.nature);
        default:
            throw std::runtime_error("Unsupported CNA type");
    }
}

bool GenomeMutations::apply(const SID& mutation, const AlleleId& allele_id)
{
    auto chr_it = find_chromosome(chromosomes, mutation.chr_id);

    return chr_it->second.apply(mutation, allele_id);
}

bool GenomeMutations::apply(const MutationSpec<SID>& mutation_spec)
{
    auto chr_it = find_chromosome(chromosomes, mutation_spec.chr_id);

    return chr_it->second.apply(mutation_spec);
}

bool GenomeMutations::remove_mutation(const GenomicPosition& genomic_position)
{
    auto chr_it = find_chromosome(chromosomes, genomic_position.chr_id);

    return chr_it->second.remove_mutation(genomic_position);
}

bool GenomeMutations::has_context_free(const SID& mutation) const
{
    auto chr_it = find_chromosome(chromosomes, mutation.chr_id);

    return chr_it->second.has_context_free(mutation);
}

bool GenomeMutations::includes(const SID& mutation) const
{
    auto chr_it = chromosomes.find(mutation.chr_id);

    if (chr_it == chromosomes.end()) {
        return false;
    }

    return chr_it->second.includes(mutation);
}

GenomeMutations GenomeMutations::copy_structure() const
{
    GenomeMutations copy;

    for (const auto& [chr_id, chromosome] : chromosomes) {
        copy.chromosomes.emplace(chr_id, chromosome.copy_structure());
    }

    return copy;
}

void GenomeMutations::duplicate_alleles()
{
    for (auto& [chr_id, chromosome] : chromosomes) {
        chromosome.duplicate_alleles();
    }
}

std::map<ChromosomeId, std::map<ChrPosition, AllelicType>>
GenomeMutations::get_allelic_types(const std::map<ChromosomeId, std::set<ChrPosition>>& break_points,
                                 const size_t& min_allelic_size) const
{
    std::map<ChromosomeId, std::map<ChrPosition, AllelicType>> allelic_map;
    for (const auto& [chr_id, chr] : get_chromosomes()) {
        allelic_map[chr_id] = chr.get_allelic_types(break_points.at(chr_id), min_allelic_size);
    }

    return allelic_map;
}

std::map<ChromosomeId, std::set<ChrPosition>>
GenomeMutations::get_CNA_break_points() const
{
    std::map<ChromosomeId, std::set<ChrPosition>> b_points;

    for (const auto& [chr_id, chr] : get_chromosomes()) {
        auto& chr_b_points = b_points[chr_id];
        for (const auto& [allele_id, allele] : chr.get_alleles()) {
            for (const auto& [pos, fragment] : allele.get_fragments()) {
                chr_b_points.insert(fragment.get_final_position());
                chr_b_points.insert(fragment.get_final_position()+1);
            }
        }
    }

    return b_points;
}

CellGenomeMutations::CellGenomeMutations():
    Cell(), GenomeMutations()
{}

CellGenomeMutations::CellGenomeMutations(const GenomeMutations& germline_mutations):
    Cell(), GenomeMutations(germline_mutations)
{}

CellGenomeMutations::CellGenomeMutations(const Mutants::Cell& cell, const GenomeMutations& genome_mutations):
    Cell(cell), GenomeMutations(genome_mutations)
{}

SampleGenomeMutations::SampleGenomeMutations(const std::string& name,
                                             const std::shared_ptr<GenomeMutations>& germline_mutations):
    germline_mutations(germline_mutations), mutations(), name(name)
{}

}   // Mutations

}   // RACES

namespace std
{

std::ostream& operator<<(std::ostream& os, const RACES::Mutations::ChromosomeMutations& chromosome_mutations)
{
    os << "Chromosome " << RACES::Mutations::GenomicPosition::chrtos(chromosome_mutations.id()) << std::endl;

    for (const auto& [allele_id, allele]: chromosome_mutations.get_alleles()) {
        os << "  " << RACES::Mutations::Allele::format_id(allele_id) << ": " << allele << std::endl;
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const RACES::Mutations::GenomeMutations& genome_mutations)
{
    for (const auto& [chr_id, chr_mutations]: genome_mutations.get_chromosomes()) {
        os << chr_mutations << std::endl;
    }

    return os;
}

}   // std
