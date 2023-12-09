/**
 * @file genotype_properties.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines genotype properties
 * @version 0.3
 * @date 2023-12-09
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

#ifndef __RACES_GENOTYPE_PROPERTIES__
#define __RACES_GENOTYPE_PROPERTIES__

#include <map>
#include <string>
#include <sstream>

#include "archive.hpp"
#include "cell_event.hpp"

namespace Races 
{

namespace Clones
{

/**
 * @brief A class to represent promoter methylation/demethylation rates
 */
class EpigeneticRates
{
    double methylation;         //!< the methylation rate
    double demethylation;       //!< the demethylation rate
    
public:
    /**
     * @brief A constructor
     * 
     * @param methylation_rate is the promoter methylation rate
     * @param demethylation_rate is the promoter demethylation rate
     */
    EpigeneticRates(const double methylation_rate, const double demethylation_rate);

    /**
     * @brief A constructor
     * 
     * @param rate is the promoter methylation/demethylation rate
     */
    explicit EpigeneticRates(const double rate);

    /**
     * @brief Get the methylation rate
     * 
     * @return a constant reference to the methylation rate
     */
    inline const double& get_methylation_rate() const
    {
        return methylation;
    }

    /**
     * @brief Set the methylation rate
     * 
     * @param rate is the new methylation rate
     * @return a reference to updated object
     */
    EpigeneticRates& set_methylation_rate(const double& rate);

    /**
     * @brief Get the demethylation rate
     * 
     * @return a constant reference to the demethylation rate
     */
    inline const double& get_demethylation_rate() const
    {
        return demethylation;
    }

    /**
     * @brief Set the demethylation rate
     * 
     * @param rate is the new demethylation rate
     * @return a reference to updated object
     */
    EpigeneticRates& set_demethylation_rate(const double& rate);
};

/**
 * @brief The methylation signature type
 * 
 * The methylation signature represents the epigenetic status 
 * of all the considered promoters. A promoter is set to be 
 * `true` when is methylated and `false` when it is not.
 */
typedef std::vector<bool> MethylationSignature;

class GenotypeProperties;

namespace Evolutions
{

    class Species;

}   // Evolutions

/**
 * @brief A class to represent species properties
 * 
 * A genotype is characterized by specific genomic (potentially 
 * unknown) mutations. The same genotype may be associated to 
 * different rates depending on its methylation signature.
 * A species can be identified by a genotype and a methylation
 * signature. 
 */
class SpeciesProperties 
{
    static unsigned int counter;                    //!< Total number of species along the computation

    SpeciesId id;                                   //!< Identifier
    GenotypeId genotype_id;                         //!< Genotype identifier

    std::string name;                               //!< Species name
    MethylationSignature methylation_signature;     //!< Methylation signature

    std::map<CellEventType, double> event_rates;    //!< Event rates

    std::map<SpeciesId, double> epigenetic_rates;   //!< Epigenetic mutation rates 

    /**
     * @brief The empty constructor
     */
    SpeciesProperties();

    /**
     * @brief A constructor
     * 
     * @param genotype is the genotype of the new object
     * @param num_of_promoters is the number of methylable promoters
     */
    SpeciesProperties(const GenotypeProperties& genotype,
                       const size_t num_of_promoters);
public:
    /**
     * @brief Get the species identifier
     * 
     * @return a constant reference to the identifier
     */
    inline const SpeciesId& get_id() const
    {
        return id;
    }

    /**
     * @brief Get the genotype identifier
     * 
     * @return a constant reference to the genotype identifier
     */
    inline const GenotypeId& get_genotype_id() const
    {
        return genotype_id;
    }

    /**
     * @brief Get the genotype name
     * 
     * @return a constant reference to the genotype name
     */
    inline const std::string& get_genotype_name() const
    {
        return name;
    }

    /**
     * @brief Get the species name
     * 
     * @return the species name
     */
    std::string get_name() const;

    /**
     * @brief Get the rate of an event 
     * 
     * @param event is the event whose rate is required
     * @return if the rate of `event` has been set, then 
     *       the rate of `event`. The value 0 otherwise
     */
    double get_rate(const CellEventType& event) const;

    /**
     * @brief Set the rate of an event 
     * 
     * @param event is the event whose rate is set
     * @param rate is the new rate for the event
     * @return a constant reference to the new rate
     */
    inline const double& set_rate(const CellEventType& event, const double rate)
    {
        return (event_rates[event] = rate);
    }

    /**
     * @brief Get event rates
     * 
     * @return a constant reference to event rates
     */
    inline const std::map<CellEventType, double>& get_rates() const
    {
        return event_rates;
    }

    /**
     * @brief Set event rates
     * 
     * @param rates are the new event rates
     */
    inline void set_rates(const std::map<CellEventType, double>& rates)
    {
        event_rates = rates;
    }

    /**
     * @brief Get the epigenetic switch rates
     * 
     * @return a constant reference to a map associating the identifier of 
     *      the species reachable with an epigenetic switch and 
     *      its rate
     */
    inline const std::map<SpeciesId, double>& get_epigenetic_switch_rates() const
    {
        return epigenetic_rates;
    }

    /**
     * @brief Get the methylation signature
     * 
     * @return the methylation signature
     */
    inline const MethylationSignature& get_methylation_signature() const
    {
        return methylation_signature;
    }

    /**
     * @brief Save the species in an archive
     * 
     * @tparam ARCHIVE is the output archive type
     * @param archive is the output archive
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::Out, ARCHIVE>, bool> = true>
    void save(ARCHIVE& archive) const
    {
        archive & id 
                & genotype_id 
                & name 
                & methylation_signature 
                & event_rates
                & epigenetic_rates;
    }

    /**
     * @brief Load a species from an archive
     * 
     * @tparam ARCHIVE is the input archive type
     * @param archive is the input archive
     * @return the loaded species 
     */
    template<typename ARCHIVE, std::enable_if_t<std::is_base_of_v<Archive::Basic::In, ARCHIVE>, bool> = true>
    static SpeciesProperties load(ARCHIVE& archive)
    {
        SpeciesProperties species_properties;

        archive & species_properties.id 
                & species_properties.genotype_id 
                & species_properties.name 
                & species_properties.methylation_signature 
                & species_properties.event_rates
                & species_properties.epigenetic_rates;

        if (SpeciesProperties::counter < static_cast<unsigned int>(species_properties.id+1)) {
            SpeciesProperties::counter = species_properties.id+1;
        }

        return species_properties;
    }

    friend class GenotypeProperties;
    friend class Evolutions::Species;
};

/**
 * @brief A class representing genotype properties
 */
class GenotypeProperties 
{
    static unsigned int counter;                        //!< Total number of genotype along the computation

    GenotypeId id;                            //!< Identifier
    std::string name;                                   //!< Genotype name

    std::vector<SpeciesProperties> species;   //!< SpeciesProperties

    template<typename SIGNATURE_TYPE>
    void validate_signature(const SIGNATURE_TYPE& signature) const;

public:

    /**
     * @brief A constructor
     * 
     * This constructor builds a genotype and its species. 
     * The number of species is determined by the size of the 
     * epigenetic event rate vector. Each cell of the vector correspond to 
     * a pair methylation/demethylation events for a specific target promoter.
     * Since each target promoter can be either methylated ("+") or 
     * non-methylated ("-"), the number of species associated to
     * a genotype is exponential in the number of the promoters and, 
     * as a consequence, in the number of elements in the epigenetic event 
     * rate vector.
     * 
     * @param name is the name of the genotype
     * @param epigenetic_event_rates is the rates of the methylation/demethylation events 
     *          on each target promoters
     */
    GenotypeProperties(const std::string& name, const std::vector<EpigeneticRates>& epigenetic_event_rates);

    /**
     * @brief A constructor
     * 
     * This constructor builds a genotype consisting in one species.
     * 
     * @param name is the name of the genotype
     */
    explicit GenotypeProperties(const std::string& name);

    /**
     * @brief Get the species identifier
     * 
     * @return a constant reference to the identifier
     */
    inline const GenotypeId& get_id() const
    {
        return id;
    }

    /**
     * @brief Get the genotype name
     * 
     * @return a constant reference to the genotype name
     */
    inline const std::string& get_name() const
    {
        return name;
    }

    /**
     * @brief Get the species associated to this genotype
     * 
     * @return a constant reference to the species
     */
    inline const std::vector<SpeciesProperties>& get_species() const
    {
        return species;
    }

    /**
     * @brief Get the number of methylable promoters in the genotype
     * 
     * @return the number of methylable promoters in the genotype
     */
    size_t num_of_promoters() const;

    /**
     * @brief Get the species associated to a methylation signature
     * 
     * @param methylation_signature is the methylation signature of the aimed species
     * @return a non-constant reference to the species associated 
     *      to `methylation_signature`
     */
    SpeciesProperties& operator[](const MethylationSignature& methylation_signature);

    /**
     * @brief Get the species associated to a methylation signature
     * 
     * @param methylation_signature is the methylation signature of the aimed species
     * @return a constant reference to the species associated 
     *      to `methylation_signature`
     */
    const SpeciesProperties& operator[](const MethylationSignature& methylation_signature) const;

    /**
     * @brief Get the species associated to a methylation signature
     * 
     * @param methylation_signature is a string representing the signature of the aimed species
     * @return a non-constant reference to the species associated 
     *      to `methylation_signature`
     */
    SpeciesProperties& operator[](const std::string& methylation_signature);

    /**
     * @brief Get the species associated to a methylation signature
     * 
     * @param methylation_signature is a string representing the signature of the aimed species
     * @return a constant reference to the species associated 
     *      to `methylation_signature`
     */
    const SpeciesProperties& operator[](const std::string& methylation_signature) const;

    static size_t string_to_index(const std::string& methylation_signature);
    static std::string index_to_string(const size_t& index, const size_t num_of_promoters);

    static size_t signature_to_index(const MethylationSignature& methylation_signature);
    static MethylationSignature index_to_signature(const size_t& index, const size_t num_of_promoters);

    inline static std::string signature_to_string(const MethylationSignature& methylation_signature)
    {
        const auto signature_index = signature_to_index(methylation_signature);

        return GenotypeProperties::index_to_string(signature_index, methylation_signature.size());
    }
};

template<typename SIGNATURE_TYPE>
void GenotypeProperties::validate_signature(const SIGNATURE_TYPE& signature) const
{
    const size_t promoters = num_of_promoters();
    if (signature.size()!=promoters) {
        std::stringstream oss;

        oss << "Wrong signature size: " 
            << get_name() << " has ";

        switch (promoters) {
            case 0:
                oss << "no methylable promoters";
                break;
            case 1:
                oss << "1 methylable promoter";
                break;
            default:
                oss << promoters << " methylable promoters";
        }

        throw std::domain_error(oss.str());
    }
}

}   // Clones

}   // Races

namespace std
{

/**
 * @brief Write information about an epigenetic rates in an output stream
 * 
 * @param out is the output stream
 * @param epigenetic_rates are the epigenetic rates to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Clones::EpigeneticRates& epigenetic_rates);

/**
 * @brief Write information about a species in an output stream
 * 
 * @param out is the output stream
 * @param genotype is the species to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Clones::SpeciesProperties& genotype);

/**
 * @brief Write information about a genotype in an output stream
 * 
 * @param out is the output stream
 * @param genotype is the genotype to be streamed
 * @return a reference to the output stream
 */
std::ostream& operator<<(std::ostream& out, const Races::Clones::GenotypeProperties& genotype);

}   // std


#endif // __RACES_GENOTYPE_PROPERTIES__