#pragma once

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <iostream>
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/polygon_soup_mesh.h"

#include "lws_options.h"
#include "tpe_energy_sc.h"
#include "tpe_flow_sc.h"
#include "Eigen/SparseLU"
#include "poly_curve_network.h"

#include "scene_file.h"

namespace LWS {
class LWSApp {
public:
    static LWSApp* instance;
    static double LWSVertexEnergy(surface::VertexPositionGeometry* geom, const surface::Vertex& vert);
    static Vector3 LWSVertexGradient(surface::VertexPositionGeometry* geom, const surface::Vertex& base,
                                     const surface::Vertex& other);
    static void DisplayWireSphere(const Vector3& center, double radius, const std::string& name);
    static void DisplayPlane(const Vector3& center, const Vector3& normal, const std::string& name);
    static void DisplayCurves(const PolyCurveNetwork* curves, const std::string& name);
    static void DisplayCyclicList(const std::vector<Vector3>& positions, const std::string& name);
    static void DisplaySobolevGradients(const PolyCurveNetwork* curves,
                                        const Eigen::MatrixXd& sobolevGradients);

    void customWindow();
    void initSolver();
    void processFileOBJ(const std::string& filename);
    void processLoopFile(const std::string& filename);
    void processSceneFile(const std::string& filename);

    void VisualizeMesh(const std::string& objName);
    void AddMeshObstacle(const std::string& objName, const Vector3& center, double p, double weight);
    void AddPlaneObstacle(const Vector3& center, const Vector3& normal, double p, double weight);
    void AddSphereObstacle(const Vector3& center, double radius);
    void SubdivideCurve();
    void MeshImplicitSurface(ImplicitSurface* surface);
    void WriteImplicitSurface();

    std::string curveName;
    PolyCurveNetwork* curves = nullptr;
    Eigen::MatrixXd sobolevGradients;
    TPEFlowSolverSC* tpeSolver = nullptr;

private:
    void centerizeLoopBarycenter(PolyCurveNetwork* curves);
    void UpdateCurvePositions();
    void outputFrame();
    void outputOBJFrame();
    void writeCurves(PolyCurveNetwork* network, const std::string& positionFilename,
                     const std::string& tangentFilename);
    void benchmarkMethods();

    SceneData sceneData;
    std::unique_ptr<surface::HalfedgeMesh> mesh;
    std::unique_ptr<surface::VertexPositionGeometry> geom;
    double initialAverageLength = 0.;
    int numStuckIterations = 0;
    int subdivideLimit = 0;
    int subdivideCount = 0;
    int currentStep = 0;
    int stepLimit = 0;
    bool perfLogging = false;
    bool useBackproj = true;
    int screenshotNum = 0;
    bool writeOBJs = false;
    int objNum = 0;
};
}
