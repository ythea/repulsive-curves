#pragma once

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "flow/gradient_constraint_enum.h"
#include "extra_potentials.h"
#include "implicit_surface.h"

namespace LWS {
struct ObstacleData {
    std::string filename;
    double weight;
};

struct PlaneObstacleData {
    Vector3 center;
    Vector3 normal;
    double weight;
};

enum class PotentialType {
    Length, LengthDiff, PinAngles, Area, VectorField
};

struct PotentialData {
    PotentialType type;
    double weight;
    std::string extraInfo;
};

class SceneData {
public:
    std::string curve_filename;
    std::vector<ObstacleData> obstacles;
    std::vector<PlaneObstacleData> planes;
    std::vector<std::string> surfacesToShow;
    std::vector<ConstraintType> constraints;
    std::vector<int> pinnedVertices;
    std::vector<int> pinnedTangents;
    std::vector<int> surfaceConstrainedVertices;
    bool pinSpecialVertices = false;
    bool pinSpecialTangents = false;
    bool pinEndpointVertices = false;
    bool constrainEndpointsToSurface = false;
    bool constrainAllToSurface = false;
    int subdivideLimit = 0;
    int iterationLimit = 0;
    double tpe_alpha = 0.;
    double tpe_beta = 0.;
    double tpe_weight = 0.;
    std::vector<PotentialData> extraPotentials;
    bool useLengthScale = false;
    double edgeLengthScale = 0.;
    bool useTotalLengthScale = false;
    double totalLengthScale = 0.;
    ImplicitSurface* constraintSurface = nullptr;

    static bool useSmoothUnion;
};

template <class Container>
void splitString(const std::string& str, Container& cont, char delim = ' ') {
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    }
}

std::string getDirectoryFromPath(std::string str);

SceneData ParseSceneFile(std::string filename);
}