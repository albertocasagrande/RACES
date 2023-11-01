/**
 * @file json_config.cpp
 * @author Alberto Casagrande (alberto.casagrande@units.it)
 * @brief Implements classes and function for reading JSON configurations
 * @version 0.7
 * @date 2023-11-01
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

#include "json_config.hpp"

#include "timed_event.hpp"

namespace Races
{

std::vector<Races::Drivers::Simulation::AxisPosition>
ConfigReader::get_corner(const nlohmann::json& sampler_region_json, const std::string& corner_field_name)
{
    if (!sampler_region_json.contains(corner_field_name)) {
        throw std::runtime_error("The \"sampler region\" must contain a "
                                    "\"" + corner_field_name + "\" field");
    }

    auto corner_json = sampler_region_json[corner_field_name];

    if (!corner_json.is_array()) {
        throw std::runtime_error("The \"" + corner_field_name + "\" must be an array of natural values");
    }

    std::vector<Races::Drivers::Simulation::AxisPosition> corner;

    for (const auto& axis_value : corner_json) {
        corner.push_back(axis_value.template get<Races::Drivers::Simulation::AxisPosition>());
    }

    if (corner.size()<2 || corner.size()>3) {
        throw std::runtime_error("\""+corner_field_name+"\" must be either a 2D or 3D array");
    }

    return corner;
}

std::pair<std::string, double>
ConfigReader::get_mutation_rate(const nlohmann::json& mutation_rate_json)
{  
    if (!mutation_rate_json.is_object()) {
        throw std::runtime_error("All the elements in \"mutation rates\" "
                                    "must be objects");
    }

    return {
        get_from<std::string>("epigenetic status", mutation_rate_json, 
                                "All the elements in \"mutation rates\""),
        get_from<double>("rate", mutation_rate_json, 
                            "All the elements in \"mutation rates\"")
    };
}

Races::Passengers::CopyNumberAlteration::Type
ConfigReader::get_CNA_type(const nlohmann::json& CNA_json)
{
    auto subtype_str = get_from<std::string>("subtype", CNA_json, "All the CNAs");
    for (auto& subtype_chr: subtype_str) {
        subtype_chr = tolower(subtype_chr);
    }

    if (subtype_str == "amplification") {
        return Races::Passengers::CopyNumberAlteration::Type::AMPLIFICATION;
    }

    if (subtype_str == "deletion") {
        return Races::Passengers::CopyNumberAlteration::Type::DELETION;
    }

    throw std::runtime_error("Unknown CNA \""+
                                CNA_json["subtype"].template get<std::string>()+"\"");
}

void
ConfigReader::add_CNA(std::list<Races::Passengers::CopyNumberAlteration>& CNAs, 
                      const nlohmann::json& CNA_json)
{
    using namespace Races::Passengers;

    auto chr_str = get_from<std::string>("chromosome", CNA_json, "All the CNAs");
    auto position = get_from<ChrPosition>("position", CNA_json, "All the CNAs");
    auto length = get_from<ChrPosition>("length", CNA_json, "All the CNAs");
    auto allele = get_from<uint8_t>("allele", CNA_json, "All the CNAs");

    GenomicPosition genomic_position(GenomicPosition::stochr(chr_str), position);

    switch (get_CNA_type(CNA_json)) {
        case CopyNumberAlteration::Type::AMPLIFICATION:
            CNAs.push_back(CopyNumberAlteration::new_amplification({genomic_position, length}, 
                                                                   allele, 0));
            break;
        case CopyNumberAlteration::Type::DELETION:
            CNAs.push_back(CopyNumberAlteration::new_deletion({genomic_position, length}, allele));
            break;
        default:
            throw std::runtime_error("Unsupported CNA type");
    }
}

void
ConfigReader::add_SNV(const std::string& driver_name, std::list<Races::Passengers::SNV>& SNVs,
                      const nlohmann::json& SNV_json)
{
    using namespace Races::Passengers;

    MutationalContext context = get_from<std::string>("context", SNV_json, 
                                                        "All the SNVs");
    auto chr_str = get_from<std::string>("chromosome", SNV_json, "All the SNVs");
    auto mutated_base = get_from<std::string>("mutated base", SNV_json, "All the SNVs");
    auto position = get_from<ChrPosition>("position", SNV_json, "All the SNVs");

    GenomicPosition genomic_position(GenomicPosition::stochr(chr_str), position);

    SNVs.emplace_back(genomic_position, context, mutated_base[0], driver_name);
}

std::map<std::string, double>
ConfigReader::get_mutation_rates(const nlohmann::json& mutation_rates_json)
{
    std::map<std::string, double> mutational_rates;

    if (!mutation_rates_json.is_array()) {
        throw std::runtime_error("The \"mutation rates\" field must be an array "
                                "of driver mutation rates");
    }

    for (const auto& mutation_rate_json : mutation_rates_json) {
        mutational_rates.insert(get_mutation_rate(mutation_rate_json));
    }

    return mutational_rates;
}

void
ConfigReader::add_driver_mutation(const std::string& driver_name,
                                  std::list<Races::Passengers::SNV>& SNVs,
                                  std::list<Races::Passengers::CopyNumberAlteration>& CNAs,
                                  const nlohmann::json& driver_mutation_json)
{
    if (!driver_mutation_json.is_object()) {
        throw std::runtime_error("All the elements in \"mutations\" must be objects");
    }

    auto type = get_from<std::string>("type", driver_mutation_json,
                                        "All the elements in \"mutations\"");
    if (type=="SNV") {
        add_SNV(driver_name, SNVs, driver_mutation_json);

        return;
    }

    if (type=="CNA") {
        add_CNA(CNAs, driver_mutation_json);

        return;
    }

    throw std::runtime_error("Unsupported mutations type \""+type+"\"");
}

double
ConfigReader::get_fraction(const nlohmann::json& fraction_json)
{
    const double fraction = fraction_json.template get<double>();

    if (fraction<0) {
        throw std::runtime_error("The \"fraction\" field must be a non-negative value");
    }

    return fraction;
}

void
ConfigReader::add_driver_mutational_properties(Races::Passengers::SpeciesMutationalProperties& mutational_properties,
                                               const Races::Drivers::Simulation::Simulation& drivers_simulation,
                                               const nlohmann::json& driver_properties_json)
{
    if (!driver_properties_json.is_object()) {
        throw std::runtime_error("All the elements in \"driver properties\" "
                                    "must be objects");
    }

    auto driver_name = get_from<std::string>("name", driver_properties_json, 
                                             "All the elements in \"driver properties\"");

    if (!driver_properties_json.contains("mutation rates")) {
        throw std::runtime_error("All the elements in \"driver properties\" "
                                "must contain a \"mutation rates\" field");
    }
    auto mutation_rates = get_mutation_rates(driver_properties_json["mutation rates"]);

    using namespace Races::Passengers;
    std::list<SNV> SNVs;
    std::list<CopyNumberAlteration> CNAs;
    if (driver_properties_json.contains("mutations")) {
        collect_driver_mutations(driver_name, SNVs, CNAs, driver_properties_json["mutations"]);
    }
    
    mutational_properties.add_species(drivers_simulation, driver_name, mutation_rates, SNVs, CNAs);
}

std::map<std::string, double>
ConfigReader::get_default_mutational_coefficients(const nlohmann::json& mutational_coeff_json)
{
    if (!mutational_coeff_json.is_object()) {
        throw std::runtime_error("The \"mutational coefficients\" field must be an object");
    }

    expecting("default", mutational_coeff_json, "The \"mutational coefficients\" field");

    return get_mutational_coefficients(mutational_coeff_json["default"]);
}

Drivers::Simulation::SampleSpecification
ConfigReader::get_sample_specification(const nlohmann::json& sample_specification_json,
                                       const std::string& field_name)
{
    if (!sample_specification_json.is_object()) {
        throw std::runtime_error("The \""+field_name+"\" field must be objects");
    }

    bool preserve_tissue=false;
    if (sample_specification_json.contains("preserve tissue")) {
        preserve_tissue = sample_specification_json["preserve tissue"].template get<bool>();
    }

    auto region = ConfigReader::get_sample_region(sample_specification_json, field_name);

    if (sample_specification_json.contains("name")) {
        auto name = sample_specification_json["name"].template get<std::string>();

        return Drivers::Simulation::SampleSpecification(name, region, preserve_tissue);
    }

    return Drivers::Simulation::SampleSpecification(region, preserve_tissue);
}

Drivers::RectangleSet
ConfigReader::get_sample_region(const nlohmann::json& sample_region_json,
                                const std::string& field_name)
{
    if (!sample_region_json.is_object()) {
        throw std::runtime_error("The \""+field_name+"\" field must be objects");
    }

    ConfigReader::expecting("lower corner", sample_region_json, field_name.c_str());
    ConfigReader::expecting("upper corner", sample_region_json, field_name.c_str());

    auto lower_vector = get_corner(sample_region_json, "lower corner");
    auto upper_vector = get_corner(sample_region_json, "upper corner");

    if (lower_vector.size() != upper_vector.size()) {
        throw std::runtime_error("The \"lower corner\" and \"upper corner\" arrays must "
                                    "have the same size");
    }

    using namespace Races::Drivers::Simulation;

    if (lower_vector.size()==2) {
        return {{lower_vector[0], lower_vector[1]},
                {upper_vector[0], upper_vector[1]}};
    }

    return {{lower_vector[0], lower_vector[1], lower_vector[2]},
            {upper_vector[0], upper_vector[1], upper_vector[2]}};
}

void 
ConfigReader::collect_driver_mutations(const std::string& driver_name,
                                        std::list<Races::Passengers::SNV>& SNVs,
                                        std::list<Races::Passengers::CopyNumberAlteration>& CNAs,
                                        const nlohmann::json& driver_mutations_json)
{
    if (!driver_mutations_json.is_array()) {
        throw std::runtime_error("The \"mutations\" field must be an array "
                                    "of mutations");
    }

    for (const auto& driver_mutation_json : driver_mutations_json) {
        add_driver_mutation(driver_name, SNVs, CNAs, driver_mutation_json);
    }
}

Races::Passengers::SpeciesMutationalProperties
ConfigReader::get_mutational_properties(const Races::Drivers::Simulation::Simulation& drivers_simulation,
                                        const nlohmann::json& simulation_json)
{
    using namespace Races::Passengers;

    SpeciesMutationalProperties mutational_properties;

    expecting("driver properties", simulation_json, "The passengers simulation configuration");

    auto& mutational_properties_json = simulation_json["driver properties"];

    if (!mutational_properties_json.is_array()) {
        throw std::runtime_error("The \"driver properties\" field must be an array "
                                    "of driver mutational properties");
    }

    for (const auto& driver_properties_json : mutational_properties_json) {
        add_driver_mutational_properties(mutational_properties, drivers_simulation,
                                         driver_properties_json);
    }

    return mutational_properties;
}

Races::Drivers::Simulation::TimedEvent
get_timed_driver_mutation(const std::map<std::string, Races::Drivers::Genotype> genotypes,
                            const nlohmann::json& timed_driver_mutation_json)
{
    using namespace Races::Drivers::Simulation;

    ConfigReader::expecting("time", timed_driver_mutation_json, "Every timed driver mutation description");

    const auto time = timed_driver_mutation_json["time"].template get<Time>();

    ConfigReader::expecting("original genotype", timed_driver_mutation_json, "Every timed driver mutation description");

    const auto& orig_genotype = genotypes.at(timed_driver_mutation_json["original genotype"].template get<std::string>());

    ConfigReader::expecting("mutated genotype", timed_driver_mutation_json, "Every timed driver mutation description");

    const auto& mutated_genotype = genotypes.at(timed_driver_mutation_json["mutated genotype"].template get<std::string>());
    
    SimulationEventWrapper driver_mutation({orig_genotype, mutated_genotype});

    return {time, driver_mutation};
}

double get_rate(const nlohmann::json& rate_json)
{
    double rate = rate_json.template get<double>();

    if (rate < 0) {
        throw std::domain_error("the rate must be non-negative");
    }

    return rate;
}

Races::Drivers::CellEventType
get_cell_event_type_by_name(const std::string& event_name)
{

    for (const auto& [type, name]: Drivers::cell_event_names) {
        if (name == event_name) {
            return type;
        }
    }

    throw std::runtime_error("Unknown cell event type \""+event_name+"\"");
}

const Races::Drivers::Simulation::Species&
find_species_by_epigenetic_name(const Races::Drivers::Simulation::Simulation& drivers_simulation,
                                const std::string& epigenetic_name)
{
    for (const auto& species : drivers_simulation.tissue()) {
        if (epigenetic_name == species.get_epigenetic_name()) {
            return species;
        }
    }

    throw std::runtime_error("Unknown species \""+epigenetic_name+"\"");
}

Races::Drivers::Simulation::TimedEvent
get_timed_rate_update(const Races::Drivers::Simulation::Simulation& drivers_simulation,
                      const nlohmann::json& timed_rate_update_json)
{
    using namespace Races::Drivers::Simulation;

    const std::string descr("Every timed rate update description");

    std::vector<std::string> fields{"time", "genotype", "status", "rate name", "rate"};

    for (const auto& field: fields) { 
        ConfigReader::expecting(field, timed_rate_update_json, descr);
    }
    const auto time = timed_rate_update_json["time"].template get<Time>();

    const auto name = (timed_rate_update_json["genotype"].template get<std::string>()
                        + timed_rate_update_json["status"].template get<std::string>());

    const Species& species = find_species_by_epigenetic_name(drivers_simulation, name);
    
    const auto rate_name = timed_rate_update_json["rate name"].template get<std::string>();
    Races::Drivers::CellEventType cell_event_type = get_cell_event_type_by_name(rate_name);

    SimulationEventWrapper rate_update({species.get_id(), cell_event_type, 
                                        get_rate(timed_rate_update_json["rate"])});

    return {time, rate_update};
}

Races::Drivers::Simulation::TimedEvent
get_timed_sampling(const nlohmann::json& timed_sampling_json)
{
    using namespace Races::Drivers::Simulation;

    const std::string descr("Every timed sampling description");

    std::vector<std::string> fields{"time", "sample"};

    for (const auto& field: fields) { 
        ConfigReader::expecting(field, timed_sampling_json, descr);
    }
    const auto time = timed_sampling_json["time"].template get<Time>();

    auto sample_specification = ConfigReader::get_sample_specification(timed_sampling_json["sample"]);

    return {time, SimulationEventWrapper(sample_specification)};
}

Races::Drivers::Simulation::TimedEvent 
ConfigReader::get_timed_event(const Races::Drivers::Simulation::Simulation& drivers_simulation,
                              const std::map<std::string, Races::Drivers::Genotype> genotypes,
                              const nlohmann::json& timed_event_json)
{
    (void)drivers_simulation;

    expecting("type", timed_event_json, "The timed event description");
    
    std::string type_name =  timed_event_json["type"].template get<std::string>();

    if (type_name == "driver mutation") {
        return get_timed_driver_mutation(genotypes, timed_event_json);
    }

    if (type_name == "liveness rate update") {
        return get_timed_rate_update(drivers_simulation, timed_event_json);
    }

    if (type_name == "sampling") {
        return get_timed_sampling(timed_event_json);
    }

    throw std::domain_error("Unsupported timed event \""+type_name+"\"");
}

void ConfigReader::expecting(const std::string& field_name, const nlohmann::json& json, 
                             const std::string& context)
{
    if (!json.contains(field_name)) {
        if (context.length()>0) {
            throw std::runtime_error("Expected a missing \""+field_name+"\" field");
        } else {
            throw std::runtime_error(context+ " must contains a \""+field_name+"\" field");
        }
    }
}

uint8_t
ConfigReader::get_number_of_alleles(const nlohmann::json& simulation_json)
{
    const double value = get_from<uint8_t>("number of alleles", simulation_json,
                                           "The passengers simulation configuration");

    if (value <= 0) {
        throw std::runtime_error("The number of allele must be positive");
    }

    return value;
}

std::map<std::string, double>
ConfigReader::get_mutational_coefficients(const nlohmann::json& mutational_coefficients_json)
{
    if (!mutational_coefficients_json.is_array()) {
        throw std::runtime_error("Mutational coefficients must be an array of object");
    }

    double total=0;

    std::map<std::string, double> mutational_coefficients;
    for (const auto& mutational_coefficient_json : mutational_coefficients_json) {

        if (!mutational_coefficient_json.is_object()) {
            throw std::runtime_error("Mutational coefficient must be an object");
        }

        const std::string SBS = get_from<std::string>("SBS", mutational_coefficient_json,
                                                        "Every mutational coefficient");
        if (mutational_coefficients.count(SBS)>0) {
            throw std::runtime_error("\""+SBS+"\" already among mutational coefficients");
        }

        double fraction = get_from<double>("fraction", mutational_coefficient_json,
                                            "Every mutational coefficient");

        total += fraction;
        if (total>1) {
            std::ostringstream oss;

            oss << "The sum of the \"fraction\" fields must be 1. "
                << "The difference is " << (1-total);

            throw std::runtime_error(oss.str());
        }

        mutational_coefficients[SBS] = fraction;
    }
    
    if (std::abs(1-total)>10*std::numeric_limits<double>::epsilon()) {
        std::ostringstream oss;

        oss << "The sum of the \"fraction\" fields must be 1. "
            << "The difference is " << (1-total);

        throw std::runtime_error(oss.str());
    }

    return mutational_coefficients;
}

}   // Races
