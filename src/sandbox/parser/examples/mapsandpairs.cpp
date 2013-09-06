/*
 * mapsandpairs.cpp
 *
 *  Created on: Sep 5, 2013
 *      Author: kglass
 */

#include <iostream>
#include <utility>
#include <map>

void printKeyValue(const std::pair<int, int> & kvp)
{
	std::cout << "<" << kvp.first << ", " << kvp.second << ">";
}

void printKeyValueSet(const std::pair<std::pair<int, int>, int > & kvp)
{
	std::cout << "<" << kvp.first.first << ", " << kvp.first.second << ">";
	std::cout << "\t" << kvp.second;
}

bool print = false;
class PairCompare
{
public:

	bool operator()(const std::pair<int, int> & a,
			const std::pair<int, int> & b)
	{

		bool rtn = a.first == b.first && a.second == b.second;
		if (print) {
		std::cout << "\t=================================" << std::endl;
		std::cout << "\ta = "; printKeyValue(a); std::cout << std::endl;
		std::cout << "\tb = "; printKeyValue(b); std::cout << std::endl;
		std::cout << "\ttrue/false = " << rtn << std::endl;
		std::cout << "\t---------------------------------" << std::endl;
		}
		return !rtn;
	}
};

class test_maps
{
public:
	void printMap() {
		for (std::map<std::pair<int, int>, int, PairCompare >::iterator it = pairMap.begin();
				it != pairMap.end(); ++it) {
			std::cout << "Pair is " << it->first.first << ", " << it->first.second << std::endl;
			std::cout << "value is " << it->second << std::endl;
		}
	}

	test_maps() {
		loadSrcMap();
		printMap();
		print = true;
		searchSrcMap();
	}

	void loadSrcMap() {
		// load a few random pairs with a few random indices
		std::pair<int, int>  dataPair;

		dataPair = std::pair<int, int>(2, 3);
		pairMap.insert(std::pair<std::pair<int, int>, int>(dataPair, 1));

		dataPair = std::pair<int, int>(4, 12);
		pairMap.insert(std::pair<std::pair<int, int>, int>(dataPair, 181));

		dataPair = std::pair<int, int>(6, 1);
		pairMap.insert(std::pair<std::pair<int, int>, int>(dataPair, 34));

		dataPair = std::pair<int, int>(4, 13);
		pairMap.insert(std::pair<std::pair<int, int>, int>(dataPair, 92));
	}

	void searchSrcMap() {
		std::cout << "++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
		std::pair<int, int>  searchPair;
		std::map<std::pair<int, int>, int>::iterator it;

		searchPair  = std::pair<int, int>(2,3);
		it = pairMap.find(searchPair);
		if (it == pairMap.end()) {
			std::cout << "item not found" << std::endl;
		} else {
			std::cout << "pair = "; printKeyValueSet(*it);std::cout << std::endl;
		}
		std::cout << "----------------------------------------------" << std::endl;

		searchPair  = std::pair<int, int>(6,1);
		it = pairMap.find(searchPair);

		if (it == pairMap.end()) {
			std::cout << "item not found" << std::endl;
		} else {
			std::cout << "pair = "; printKeyValueSet(*it);std::cout << std::endl;
		}
		std::cout << "----------------------------------------------" << std::endl;

		searchPair  = std::pair<int, int>(4,13);
		it = pairMap.find(searchPair);

		if (it == pairMap.end()) {
			std::cout << "item not found" << std::endl;
		} else {
			std::cout << "pair = "; printKeyValueSet(*it);std::cout << std::endl;
		}
		std::cout << "----------------------------------------------" << std::endl;

		searchPair  = std::pair<int, int>(4,12);
		it = pairMap.find(searchPair);

		if (it == pairMap.end()) {
			std::cout << "item not found" << std::endl;
		} else {
			std::cout << "pair = "; printKeyValueSet(*it);std::cout << std::endl;
		}
		std::cout << "----------------------------------------------" << std::endl;

		searchPair  = std::pair<int, int>(12,1);
		it = pairMap.find(searchPair);

		if (it == pairMap.end()) {
			std::cout << "item not found" << std::endl;
		} else {
			std::cout << it->first.first << " " << it->first.second << std::endl;
			std::cout << it->second << " should be 1" << std::endl;
		}
		std::cout << "++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

	}

	std::map<std::pair<int, int>, int>  pairMap;
//	std::map<std::pair<int, int>, int, PairCompare >  pairMap;
};


int main()
{
	test_maps m;
};
