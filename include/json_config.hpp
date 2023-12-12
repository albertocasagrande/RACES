/**
 * @file json_config.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes and function for reading JSON configurations
 * @version 0.11
 * @date 2023-12-12
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

#ifndef __RACES_JSON_CONFIG__
#define __RACES_JSON_CONFIG__

#include <map>
#include <vector>

#include <nlohmann/json.hpp>

#include "simulation.hpp"
#include "position_set.hpp"
#include "mutation_engine.hpp"
#include "timed_event.hpp"

namespace Races
{

class ConfigReader
{
    /**
     * @brief Get a sample region corner
     *
     * @param sampler_region_json is the JSON of the region to sample
     * @param corner_field_name is the corner field name
     * @return std::vector<Races::Mutants::Evolutions::AxisPosition>
     */
    static std::vector<Races::Mutants::Evolutions::AxisPosition>
    get_corner(const nlohmann::json& sampler_region_json, const std::string& corner_field_name);

    /**
     * @brief Get a mutation rate
     *
     * @param mutation_rate_json is the JSON of the mutation rate
     * @return a pair epigenetic status/mutation rate
     */
    static std::pair<std::string, double>
    get_mutation_rate(const nlohmann::json& mutation_rate_json);

    /**
     * @brief Get the CNA type
     *
     * @param CNA_json is the JSON of the CNA
     * @return the type of the CNA
     */
    static Races::Mutations::CopyNumberAlteration::Type
    get_CNA_type(const nlohmann::json& CNA_json);

    /**
     * @brief Add a CNA to a list
     *
     * @param CNAs is the list of CNA in which the read CNA should be inserted
     * @param CNA_json is the JSON of the CNA to be inserted into the list
     */
    static void
    add_CNA(std::list<Races::Mutations::CopyNumberAlteration>& CNAs,
            const nlohmann::json& CNA_json);

    /**
     * @brief Add a SNV to a list
     *
     * @param mutant_name is the mutant name
     * @param SNVs is the list of SNV in which the read SNV should be inserted
     * @param SNV_json is the JSON of the SNV to be inserted into the list
     */
    static void
    add_SNV(const std::string& mutant_name, std::list<Races::Mutations::SNV>& SNVs,
            const nlohmann::json& SNV_json);

    /**
     * @brief Get the mutation rates
     *
     * @param mutation_rates_json is the JSON of the mutation rates
     * @return a map associating an epigenetic status to its mutation rate
     */
    static std::map<std::string, double>
    get_mutation_rates(const nlohmann::json& mutation_rates_json);

    /**
     * @brief Add the mutations to SNV and CNA lists
     *
     * @param mutant_name is the mutant name
     * @param SNVs is the list of SNVs in which the new SNVs must be inserted
     * @param CNAs is the list of CNAs in which the new CNAs must be inserted
     * @param mutation_json is the JSON of the mutation
     */
    static void
    schedule_mutation(const std::string& mutant_name,
                      std::list<Races::Mutations::SNV>& SNVs,
                      std::list<Races::Mutations::CopyNumberAlteration>& CNAs,
                      const nlohmann::json& mutation_json);

    /**
     * @brief Get the fraction field
     *
     * @param fraction_json is the JSON of the faction field
     * @return the fraction value
     */
    static double get_fraction(const nlohmann::json& fraction_json);

    /**
     * @brief Add mutational properties
     *
     * @param mutational_properties are the mutational properties of all the species
     * @param simulation is the simulation
     * @param mutational_properties_json is the JSON of the properties
     */
    static void
    add_mutational_properties(Races::Mutations::MutationalProperties& mutational_properties,
                              const Races::Mutants::Evolutions::Simulation& simulation,
                              const nlohmann::json& mutational_properties_json);

public:
    /**
     * @brief Get the default mutational coefficients
     *
     * @param mutational_coeff_json is the JSON of the mutational coefficients
     * @return a map associating to a set of mutational signatures their percentage
     *          in the default mutational configuration
     */
    static std::map<std::string, double>
    get_default_mutational_coefficients(const nlohmann::json& mutational_coeff_json);

    /**
     * @brief Add the timed mutational coefficients to a mutation engine
     *
     * @tparam GENOME_WIDE_POSITION is the genome wide position type
     * @tparam RANDOM_GENERATOR is the type of the random generator
     * @param engine is a mutation engine
     * @param mutational_coeff_json is the JSON of the mutational coefficients
     */
    template<typename GENOME_WIDE_POSITION, typename RANDOM_GENERATOR>
    static void
    add_timed_mutational_coefficients(Races::Mutations::MutationEngine<GENOME_WIDE_POSITION,RANDOM_GENERATOR>& engine,
                                      const nlohmann::json& mutational_coeff_json)
    {
        std::map<double, std::map<std::string, double>> timed_coefficients;

        if (!mutational_coeff_json.is_object()) {
            throw std::runtime_error("The \"mutational coefficients\" field must be an object");
        }

        if (mutational_coeff_json.contains("timed")) {

            auto& timed_coeff_json = mutational_coeff_json["timed"];

            if (!timed_coeff_json.is_array()) {
                throw std::runtime_error("The optional \"timed\" field must be an array");
            }

            for (const auto& timed_json : timed_coeff_json) {
                if (!timed_json.is_object()) {
                    throw std::runtime_error("The elements of the \"timed\" field must be objects");
                }

                double time = get_from<double>("time", timed_json,
                                               "The elements of the \"timed\" field");

                if (timed_coefficients.count(time)>0) {
                    throw std::runtime_error("Two elements of the \"timed\" field have the "
                                             "same \"time\"");
                }

                expecting("coefficients", timed_json, "The elements of the \"timed\" field");

                engine.add(time, get_mutational_coefficients(timed_json["coefficients"]));
            }
        }
    }

    /**
     * @brief Get a sample region
     *
     * @param sample_region_json is the sample region JSON
     * @param field_name is the name of the field
     * @return the rectangle position set
     */
    static Mutants::RectangleSet
    get_sample_region(const nlohmann::json& sample_region_json,
                      const std::string& field_name="sample region");

    /**
     * @brief Get a sample specification
     *
     * @param sample_specification_json is the sample specification JSON
     * @param field_name is the name of the field
     * @return a sample specification
     */
    static Mutants::Evolutions::SampleSpecification
    get_sample_specification(const nlohmann::json& sample_specification_json,
                             const std::string& field_name="sample");

    /**
     * @brief Collect the mutations
     *
     * @param[in] mutant_name is the mutant name
     * @param[in,out] SNVs is the list of mutant SNVs
     * @param[in,out] CNAs is the list of mutant CNAs
     * @param[in] mutations_json is the JSON of mutations
     */
    static void
    collect_mutations(const std::string& mutant_name,
                      std::list<Races::Mutations::SNV>& SNVs,
                      std::list<Races::Mutations::CopyNumberAlteration>& CNAs,
                      const nlohmann::json& mutations_json);

    /**
     * @brief Get the species mutational properties
     *
     * @param simulation is a simulation
     * @param configuration_json is the JSON of the simulation configuration
     * @return the species mutational properties
     */
    static Races::Mutations::MutationalProperties
    get_mutational_properties(const Races::Mutants::Evolutions::Simulation& simulation,
                              const nlohmann::json& configuration_json);

    /**
     * @brief Extract a timed event from a JSON object
     *
     * @param simulation is a simulation
     * @param name2mutant is the map from name to mutant
     * @param timed_event_json is the JSON of the timed event
     * @return the time event described in `timed_event_json`
     */
    static Races::Mutants::Evolutions::TimedEvent
    get_timed_event(const Races::Mutants::Evolutions::Simulation& simulation,
                    const std::map<std::string, Races::Mutants::MutantProperties> name2mutant,
                    const nlohmann::json& timed_event_json);

    /**
     * @brief Raise an exception if a field is not available in a JSON
     *
     * @param field_name is the aimed field name
     * @param json is the JSON that must contain the field
     * @param context is a description of the context
     */
    static
    void expecting(const std::string& field_name, const nlohmann::json& json,
                   const std::string& context="");

    /**
     * @brief Get a value of a JSON field
     *
     * @tparam TYPE is the type of the JSON field value
     * @param field_name is the aimed field name
     * @param json is the JSON that must contain the field
     * @param context is a description of the context
     * @return the value of the JSON field `field_name`
     */
    template<typename TYPE>
    static TYPE get_from(const std::string& field_name, const nlohmann::json& json,
                         const std::string& context="")
    {
        expecting(field_name, json, context);

        return json[field_name].template get<TYPE>();
    }

    /**
     * @brief Get the number of wild-type genome alleles from the configuration
     *
     * @param simulation_json is the JSON simulation configuration
     * @param chromosome_ids is the vector of the genome chromosome ids
     * @return a pair whose first element is the default number of allele per
     *          chromosome and the second element is a map reporting the exceptions
     *          The map associates the chromosome names to the corresponding
     *          number of alleles
     */
    static std::map<Mutations::ChromosomeId, size_t>
    get_number_of_alleles(const nlohmann::json& simulation_json,
                          const std::vector<Mutations::ChromosomeId>& chromosome_ids);

    /**
     * @brief Get the mutational coefficients
     *
     * @param mutational_coefficients_json the JSON of the mutational coefficients
     * @return a map associating a set of mutational signature to their percentage
     *          in the mutational coefficients
     */
    static std::map<std::string, double>
    get_mutational_coefficients(const nlohmann::json& mutational_coefficients_json);
};

}   // Races

#endif // __RACES_JSON_CONFIG__
