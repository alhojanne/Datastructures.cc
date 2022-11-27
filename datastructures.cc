#include "datastructures.hh"
#include <memory>
#include <random>
#include <cmath>
#include <algorithm>


/**
 * @brief Datastructures::Datastructures initializes Datastructures-class
 */
Datastructures::Datastructures()
{
}

/**
 * @brief Datastructures::~Datastructures Destructor of Datastructures class
 */
Datastructures::~Datastructures()
{
}

/**
 * @brief Datastructures::station_count Calculates how many stations are in the datastructure
 * @return number of stations
 */
unsigned int Datastructures::station_count()
{
    return stations.size();
}

/**
 * @brief Datastructures::clear_all removes all stations and regions from the datastructures used
 */
void Datastructures::clear_all()
{
    stations.clear();
    stations_vector.clear();
    regions.clear();
}

/**
 * @brief Datastructures::all_stations method lists id's of all stations
 * @return a vector of all stations' id's
 */
std::vector<StationID> Datastructures::all_stations()
{
    std::vector<StationID> station_ids;
    for (auto& iter: stations)
    {
        station_ids.push_back(iter.first);
    }
    return station_ids;
}

/**
 * @brief Datastructures::add_station Adds a new station
 * @param id Added stations' id
 * @param name Added stations' name
 * @param xy Added station's coordinates
 * @return true if station added, false if station already existed
 */
bool Datastructures::add_station(StationID id, const Name& name, Coord xy)
{
    if (stations.find(id) != stations.end())
    {
        return false;
    }
    Station added_station = Station(id, name, xy);
    stations.insert(std::make_pair(id, added_station));
    stations_vector.push_back(added_station);
    //stations_alphabetical.insert(std::make_pair(name, added_station));
    return true;
}

/**
 * @brief Datastructures::get_station_name finds stations' name by given id
 * @param id id of the station to find
 * @return station's name or NO_NAME if station not found
 */
Name Datastructures::get_station_name(StationID id)
{
    auto found_station = find_element(id);
    if (found_station != stations.end())
    {
        return stations.at(id).station_name;
    }
    return NO_NAME;
}

/**
 * @brief Datastructures::get_station_coordinates finds station's coordinates by given id
 * @param id id of the station to find
 * @return found stations' coordinates or NO_COORD if station not found
 */
Coord Datastructures::get_station_coordinates(StationID id)
{
    auto found_station = find_element(id);
    if (found_station != stations.end())
    {
        return stations.at(id).coord;
    }
    return NO_COORD;
}

/**
 * @brief Datastructures::stations_alphabetically Lists all stations in alphabetical order
 * @return a vector of all station id's sorted in alphabetical order
 */
std::vector<StationID> Datastructures::stations_alphabetically()
{
    std::vector<StationID> station_ids;
    auto compare_name = [](auto x, auto y)
    {
        return x.station_name < y.station_name;
    };
    if (!std::is_sorted(stations_vector.begin(), stations_vector.end(), compare_name))
    {
        std::sort(stations_vector.begin(), stations_vector.end(), compare_name);
    }
    for (auto& iter : stations_vector)
    {
        station_ids.push_back(iter.id);
    }
    return station_ids;
}

/**
 * @brief Datastructures::stations_distance_increasing Lists all stations ordered by their distance from the origin (0,0)
 * @return a vector of all stations id sorted by distance
 */
std::vector<StationID> Datastructures::stations_distance_increasing()
{
    auto compare_distance_coords = [](auto x, auto y)
    {
        if (x.distance == y.distance)
        {
            return x.coord.y < y.coord.y;
        }
        else
        {
            return x.distance < y.distance;
        }
    };
    if (!std::is_sorted(stations_vector.begin(), stations_vector.end(), compare_distance_coords))
    {
        std::sort(stations_vector.begin(), stations_vector.end(), compare_distance_coords);
    }
    std::vector<StationID> station_ids_distance_increasing;
    station_ids_distance_increasing.reserve(stations.size());
    for (auto& iter : stations_vector)
    {
        station_ids_distance_increasing.push_back(iter.id);
    }
    return station_ids_distance_increasing;

}

/**
 * @brief Datastructures::find_station_with_coord Finds a station by a given coordinate
 * @param xy coordinate of the station
 * @return found station's id or NO_STATION if no station found
 */
StationID Datastructures::find_station_with_coord(Coord xy)
{
    for (auto& iter : stations_vector)
    {
        if (iter.coord == xy)
        {
            return iter.id;
        }
    }
    return NO_STATION;
}

/**
 * @brief Datastructures::change_station_coord Changes station's coordinates
 * @param id station's id
 * @param newcoord station's new coordinates
 * @return true if coordinates changed, false if no station found
 */
bool Datastructures::change_station_coord(StationID id, Coord newcoord)
{
    auto found_station = find_element(id);
    auto found_station_in_vector = std::find_if(stations_vector.begin(), stations_vector.end(), [id](auto x){ return x.id == id;
    });
    if (found_station != stations.end())
    {
        found_station->second.coord = newcoord;
        found_station_in_vector->coord = newcoord;
        found_station_in_vector->distance = std::sqrt(std::pow(found_station_in_vector->coord.x, 2) + std::pow(found_station_in_vector->coord.y, 2));
        return true;
    }
    return false;
}

/**
 * @brief Datastructures::add_departure Adds a departure of a train to a station
 * @param stationid station to add the departure to
 * @param trainid departing train's id
 * @param time time of the departure
 * @return false if no station found, true if departure added
 */
bool Datastructures::add_departure(StationID stationid, TrainID trainid, Time time)
{
    auto found_station = find_element(stationid);

    if (found_station == stations.end())
    {
        return false;
    }
    else if (stations.at(stationid).trains.find(time) == stations.at(stationid).trains.end())
    {
        std::vector<TrainID> added_train;
        added_train.push_back(trainid);
        stations.at(stationid).trains.insert(std::make_pair(time, added_train));
        return true;
    }
    for (auto& iter : stations.at(stationid).trains.at(time))
    {
        if (iter == trainid)
        {
            return false;
        }
    }
    stations.at(stationid).trains.at(time).push_back(trainid);
    return true;

}

/**
 * @brief Datastructures::remove_departure Removes a departure of a train from a station
 * @param stationid station's id
 * @param trainid removed train's id
 * @param time time of the removed departure
 * @return true if departure removed, false if no mathchin station or train found
 */
bool Datastructures::remove_departure(StationID stationid, TrainID trainid, Time time)
{
    auto found_station = find_element(stationid);
    if (found_station == stations.end())
    {
        return false;
    }
    auto found_departure_time = (stations.at(stationid).trains.find(time));
    if (found_departure_time == stations.at(stationid).trains.end())
    {
        return false;
    }
    std::vector<TrainID> possible_trains = stations.at(stationid).trains.at(time);
    auto found_train = std::find(possible_trains.begin(), possible_trains.end(), trainid);
    if (found_train == possible_trains.end())
    {
        return false;
    }
    else
    {
        possible_trains.erase(found_train);
        stations.at(stationid).trains.at(time) = possible_trains;
        return true;
    }
    return false;
}

/**
 * @brief Datastructures::station_departures_after Lists all departures from a station after given time
 * @param stationid station's id
 * @param time time of day
 * @return  {NOTIME, NO_TRAIN} if no station found, list of all departures after
 * given time from the station otherwise
 */
std::vector<std::pair<Time, TrainID>> Datastructures::station_departures_after(StationID stationid, Time time)
{
    auto found_station = find_element(stationid);
    std::vector<std::pair<Time, TrainID>> departures;
    if (found_station == stations.end())
    {
        departures.push_back({NO_TIME, NO_TRAIN});
        return departures;
    }

    auto first_time = stations.at(stationid).trains.find(time);
    if (first_time == stations.at(stationid).trains.end())
    {
        first_time = stations.at(stationid).trains.upper_bound(time);
    }
    auto extract_trains = [&departures](auto x)
    {
        std::sort(x.second.begin(), x.second.end());
        for (auto& iter : x.second)
        {
            departures.push_back(std::make_pair(x.first, iter));
        }
    };
    std::for_each(first_time, stations.at(stationid).trains.end(), extract_trains);
    return departures;
}

/**
 * @brief Datastructures::add_region Adds a new region
 * @param id region's id
 * @param name region's name
 * @param coords region's coordinates
 * @return false if region already existed, true after adding the new region
 */
bool Datastructures::add_region(RegionID id, const Name &name, std::vector<Coord> coords)
{
    auto found_region = find_region(id);
    if (found_region != regions.end())
    {
        return false;
    }
    Region added_region = Region(id, name, coords);
    regions.insert(std::make_pair(id, added_region));
    return true;
}

/**
 * @brief Datastructures::all_regions Lists all region's ids
 * @return a vector of ids of all regions'
 */
std::vector<RegionID> Datastructures::all_regions()
{
    std::vector<RegionID> region_ids;
    for (auto& iter: regions)
    {
        region_ids.push_back(iter.first);
    }
    return region_ids;
}

/**
 * @brief Datastructures::get_region_name Finds a region by id
 * @param id region's id
 * @return NO_NAME if no region found, region's name if found
 */
Name Datastructures::get_region_name(RegionID id)
{
    auto found_region = find_region(id);
    if (found_region == regions.end())
    {
        return NO_NAME;
    }
    return regions.at(id).region_name;
}

/**
 * @brief Datastructures::get_region_coords Finds region's coordinates by it's id
 * @param id region's id
 * @return NO_COORD if no region found, region's coordinates if found
 */
std::vector<Coord> Datastructures::get_region_coords(RegionID id)
{
    auto found_region = find_region(id);
    std::vector<Coord> no_region;
    if (found_region == regions.end())
    {
        no_region.push_back(NO_COORD);
        return no_region;
    }
    return regions.at(id).region_coords;
}

/**
 * @brief Datastructures::add_subregion_to_region adds a subregion to a region
 * @param id subregion's id
 * @param parentid parent region's id
 * @return true if subregion added to region, false if region or subregion not found
 */
bool Datastructures::add_subregion_to_region(RegionID id, RegionID parentid)
{
    auto found_region = find_region(id);
    auto found_sub_region = find_region(parentid);
    if (found_region == regions.end() || found_sub_region == regions.end())
    {
        return false;
    }
    else if (regions.at(id).parent_region != nullptr)
    {
        return false;
    }
    else
    {
        regions.at(id).parent_region = &regions.at(parentid);
        regions.at(parentid).sub_regions.push_back(&regions.at(id));
        return true;
    }
}

/**
 * @brief Datastructures::add_station_to_region Adds a station to a region
 * @param id id of the station
 * @param parentid id of the region
 * @return true if station added to region, false if station or region not found or
 * station already belongs to another region
 */
bool Datastructures::add_station_to_region(StationID id, RegionID parentid)
{
    auto found_region = find_region(parentid);
    auto found_station = find_element(id);
    if (found_region == regions.end() || found_station == stations.end())
    {
        return false;
    }
    else if (found_station->second.stations_region != nullptr)
    {
        return false;
    }
    else
    {
        found_station->second.stations_region = &found_region->second;
        return true;
    }
}

/**
 * @brief Datastructures::station_in_regions finds the region of the station and
 * all parent regions the region belongs to
 * @param id station's id
 * @return a vector of all the found regions or NO_REGION if
 * station not found or station doesn't belong to any region
 */
std::vector<RegionID> Datastructures::station_in_regions(StationID id)
{
    std::vector<RegionID> regions;
    auto found_station = find_element(id);
    auto has_region = found_station->second.stations_region;
    if (found_station == stations.end() || has_region == nullptr)
    {
        regions.push_back(NO_REGION);
        return regions;
    }
    found_region(found_station->second.stations_region->region_id, regions);
    return regions;
}

/**
 * @brief Datastructures::all_subregions_of_region finds all subregion of a region
 * @param id id of the region
 * @return a vector of all subregions's ids of the given region or
 * NO_REGION if no region found
 */
std::vector<RegionID> Datastructures::all_subregions_of_region(RegionID id)
{
    std::vector<RegionID> regions_vector;
    auto found_region = find_region(id);;
    if (found_region == regions.end())
    {
        regions_vector.push_back(NO_REGION);
        return regions_vector;
    }
    find_child_regions(id, regions_vector);
    return regions_vector;

}

/**
 * @brief Datastructures::stations_closest_to finds up to 3 closest stations to a given coordinate
 * @param xy given coordinate
 * @return a vector of stations' ids closest to the given coordinate
 */
std::vector<StationID> Datastructures::stations_closest_to(Coord xy)
{
    std::vector<std::pair<StationID, Distance>> distances_from_station;
    distances_from_station.reserve(stations_vector.size());
    for (auto& iter : stations_vector)
    {
        Distance distance = ((xy.x - iter.coord.x) * (xy.x - iter.coord.x) + (xy.y - iter.coord.y) *(xy.y - iter.coord.y) );

        auto distance_id_pair = std::make_pair(iter.id, distance);
        distances_from_station.push_back(distance_id_pair);
    }
    auto iter = distances_from_station.begin();
    std::vector<std::pair<StationID, Distance>> found_stations;

    while (iter < distances_from_station.end() && found_stations.size() < 3)
    {
        auto min = std::min_element(iter, distances_from_station.end(), [](auto x, auto y){return x.second < y.second;});
        found_stations.push_back(*min);
        //std::move(min, min,distances_from_station.begin());
        distances_from_station.erase(min);

    }
    std::vector<StationID> results;
    for (auto& iter : found_stations)
    {
        results.push_back(iter.first);
    }
    return results;
}

/**
 * @brief Datastructures::remove_station removes a station
 * @param id id of the station
 * @return true if station removed, false if no station found
 */
bool Datastructures::remove_station(StationID id)
{
    if (find_element(id) == stations.end())
    {
        return false;
    }
    stations.erase(id);
    auto found_station = std::find_if(stations_vector.begin(), stations_vector.end(), [id](auto x){return x.id == id;});
    stations_vector.erase(found_station);
    return true;
}

/**
 * @brief Datastructures::common_parent_of_regions find a closest common parent region
 * of 2 regions
 * @param id1 first region's id
 * @param id2 second region's id
 * @return  closest common parent region's id or NO_REGION
 */
RegionID Datastructures::common_parent_of_regions(RegionID id1, RegionID id2)
{
    auto found_first = find_region(id1);
    auto found_second = find_region(id2);
    std::vector<RegionID> first_regions_parents;
    std::vector<RegionID> second_regions_parents;
    auto first_regions_parent = found_first->second.parent_region;
    auto second_regions_parent = found_second->second.parent_region;
    if (found_first == regions.end() || found_second == regions.end())
    {
        return NO_REGION;
    }
    if (first_regions_parent == nullptr || second_regions_parent == nullptr)
    {
        return NO_REGION;
    }
    found_region(found_first->second.parent_region->region_id, first_regions_parents);
    found_region(found_second->second.parent_region->region_id, second_regions_parents);
    std::vector<RegionID> result_vector;
    auto found_common_region = std::find_first_of(first_regions_parents.begin(), first_regions_parents.end(), second_regions_parents.begin(), second_regions_parents.end());
    if (found_common_region == first_regions_parents.end())
    {
        return NO_REGION;
    }
    return *found_common_region;
}

/**
 * @brief Datastructures::find_element finds a given station
 * @param id station's id
 * @return iterator to the found station
 */
Datastructures::mapIterator Datastructures::find_element(StationID id)
{
    return stations.find(id);
}

/**
 * @brief Datastructures::find_region finds a given region
 * @param id region's id
 * @return iterator to the found region
 */
Datastructures::regionIterator Datastructures::find_region(RegionID id)
{
    return regions.find(id);
}

/**
 * @brief Datastructures::found_region finds region's all parent regions
 * and pushes them to a vector
 * @param parentid first parent region's id
 * @param regions_in_region a vector of all the found regions
 */
void Datastructures::found_region(RegionID parentid, std::vector<RegionID>& regions_in_region)
{
    auto next_region = regions.at(parentid);
    regions_in_region.push_back(parentid);
    if (next_region.parent_region != nullptr)
    {
        next_region = *regions.at(parentid).parent_region;
        found_region(next_region.region_id, regions_in_region);
    }
}

/**
 * @brief Datastructures::find_child_regions finds all child regions of a given region
 * @param id id of the region
 * @param child_regions empty vector used to store the child regions
 */
void Datastructures::find_child_regions(RegionID id, std::vector<RegionID> &child_regions)
{
    for (auto& iter : regions.at(id).sub_regions)
    {
        child_regions.push_back(iter->region_id);
        find_child_regions(iter->region_id, child_regions);
    }
}



