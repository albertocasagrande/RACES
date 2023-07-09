/**
 * @file driver.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Driver genotype representation
 * @version 0.4
 * @date 2023-07-09
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

#include <iostream>
#include <sstream>

#include "driver_genotype.hpp"
#include "cell_event.hpp"

namespace Races {

unsigned int Races::EpigeneticGenotype::counter = 0;
unsigned int Races::SomaticGenotype::counter = 0;

EpigeneticRates::EpigeneticRates(const double methylation_rate, const double demethylation_rate):
        std::pair<double,double>(methylation_rate, demethylation_rate)
{}

EpigeneticRates::EpigeneticRates(const double rate):
        EpigeneticRates(rate, rate)
{}

std::ostream& operator<<(std::ostream& out, const EpigeneticRates& epigentic_rates)
{
    out << "{methylation: " << epigentic_rates.methylation_rate() 
            << ",demethylation: " << epigentic_rates.demethylation_rate() << "}";
    return out;
}

EpigeneticGenotype::EpigeneticGenotype():
    id(0), somatic_id(0)
{}

EpigeneticGenotype::EpigeneticGenotype(const SomaticGenotype& somatic_genotype,
                                       const size_t num_of_promoters):
    id(counter++), somatic_id(somatic_genotype.get_id())
{
    const size_t epigenetic_index = somatic_genotype.epigenetic_genotypes().size();

    name = somatic_genotype.get_name() + SomaticGenotype::index_to_string(epigenetic_index, num_of_promoters);

    methylation_signature = SomaticGenotype::index_to_signature(epigenetic_index, num_of_promoters);
}

std::ostream& operator<<(std::ostream& out, const EpigeneticGenotype& genotype)
{
    out << "{name: \""<< genotype.get_name() << "\", id: " << genotype.get_id() 
        << ", event_rates: {";

    std::string sep="";
    for (const auto& [event, rate]: genotype.get_rates()) {
        out << sep << cell_event_names[event] << ": " << rate;
        sep = ", ";
    }

    const auto& e_rates = genotype.get_epigenentic_mutation_rates();
    if (e_rates.size()>0) {
        out << "}, epigenetic rates: {";
        sep="";
        for (const auto& [dst_id, rate]: e_rates) {
            out << sep << dst_id << ": " << rate;
            sep = ", ";
        }
        out << "}";
    }

    out << "}";

    return out;
}

SomaticGenotype::SomaticGenotype(const std::string& name,
                               const std::vector<EpigeneticRates> epigenetic_event_rates):
    id(counter++), name(name)
{
    size_t epigenetic_mutations = 1<<epigenetic_event_rates.size();

    for (size_t i=0; i<epigenetic_mutations; ++i) {
        e_genotypes.push_back(EpigeneticGenotype(*this, epigenetic_event_rates.size()));
    }

    for (auto& e_genotype: e_genotypes) {
        auto e_signature = e_genotype.get_methylation_signature();
        auto& e_rates = e_genotype.epigenetic_rates;
        
        for (size_t i=0; i<epigenetic_event_rates.size(); ++i) {
            // get the signature of the genotype reachable by 
            // methylating/demethylating the i-th promoter
            e_signature[i] = !e_signature[i];

            // get the index of the genotype reachable by 
            // methylating/demethylating the i-th promoter
            size_t index = signature_to_index(e_signature);

            // get the identifier of the genotype reachable by 
            // methylating/demethylating the i-th promoter
            EpigeneticGenotypeId dst_id = e_genotypes[index].get_id();

            if (e_signature[i]) {  // if the promoter is methylated
                // a methylation event must occur
                e_rates[dst_id] = epigenetic_event_rates[i].methylation_rate();
            } else {
                // a demethylation event must occur
                e_rates[dst_id] = epigenetic_event_rates[i].demethylation_rate();
            }

            // revert to the original signature
            e_signature[i] = !e_signature[i];
        }
    }
}

SomaticGenotype::SomaticGenotype(const std::string& name):
    SomaticGenotype(name, std::vector<EpigeneticRates>())
{}

size_t SomaticGenotype::num_of_promoters() const
{
    if (e_genotypes.size()==0) {
        return 0;
    }

    return e_genotypes[0].get_methylation_signature().size();
}

EpigeneticGenotype& SomaticGenotype::operator[](const MethylationSignature& methylation_signature)
{
    validate_signature(methylation_signature);

    return e_genotypes[SomaticGenotype::signature_to_index(methylation_signature)];
}

const EpigeneticGenotype& SomaticGenotype::operator[](const MethylationSignature& methylation_signature) const
{
    validate_signature(methylation_signature);

    return e_genotypes[SomaticGenotype::signature_to_index(methylation_signature)];
}

EpigeneticGenotype& SomaticGenotype::operator[](const std::string& methylation_signature)
{
    validate_signature(methylation_signature);

    return e_genotypes[SomaticGenotype::string_to_index(methylation_signature)];
}

const EpigeneticGenotype& SomaticGenotype::operator[](const std::string& methylation_signature) const
{
    validate_signature(methylation_signature);

    return e_genotypes[SomaticGenotype::string_to_index(methylation_signature)];
}

size_t SomaticGenotype::string_to_index(const std::string& methylation_signature)
{
    size_t index=0, digit_value=(1 << methylation_signature.size());

    for (const char& c: methylation_signature) {
        digit_value >>= 1;
        if (c=='+') {
            index |= digit_value;
        }
    }

    return index;
}

std::string SomaticGenotype::index_to_string(const size_t& index, const size_t num_of_promoters)
{
    size_t digit_value=(1 << num_of_promoters);
    std::ostringstream oss;

    while (digit_value>1) {
        digit_value >>= 1;
        oss << (((digit_value&index) != 0) ? "+" : "-");
    }

    return oss.str();
}

size_t SomaticGenotype::signature_to_index(const MethylationSignature& methylation_signature)
{
    size_t index=0, digit_value=(1 << methylation_signature.size());

    for (const bool& b_value: methylation_signature) {
        digit_value >>= 1;
        if (b_value) {
            index |= digit_value;
        }
    }

    return index;
}

MethylationSignature SomaticGenotype::index_to_signature(const size_t& index, const size_t num_of_promoters)
{
    MethylationSignature signature;
    size_t digit_value=(1 << num_of_promoters);

    while (digit_value>1) {
        digit_value >>= 1;
        signature.push_back((digit_value&index) != 0);
    }

    return signature;
}

std::ostream& operator<<(std::ostream& out, const SomaticGenotype& genotype)
{
    out << "{name: \"" << genotype.get_name() << "\", id: " 
        << genotype.get_id() << ", epigenetic genotypes=[";

    const auto& e_genotypes = genotype.epigenetic_genotypes();
    
    std::string sep = "";
    if (e_genotypes.size()>1) {
        sep = "\n";
    }

    for (const auto& e_genotype: e_genotypes) {
        out << sep << e_genotype;
        sep = ",\n";
    }

    out << "]}";
    return out; 
}

};