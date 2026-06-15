#pragma once

#include <vector>
#include "geometrycentral/utilities/vector3.h"
#include "Eigen/Core"
#include "utils.h"
#include "libgmultigrid/multigrid_operator.h"
#include <unordered_set>
#include "flow/gradient_constraint_enum.h"
#include "implicit_surface.h"
#include "multigrid/constraint_projector_operator.h"

namespace LWS {
using namespace geometrycentral;
class PolyCurveNetwork;
struct CurveVertex;

struct CurveEdge {
    CurveVertex* nextVert = nullptr;
    CurveVertex* prevVert = nullptr;
    int id = 0;
    int component = 0;

    inline CurveVertex* Opposite(const CurveVertex* v) const;
    inline Vector3 Vector() const;
    inline Vector3 Tangent() const;
    inline double Length() const;
    inline Vector3 Midpoint() const;
    //inline Vector3 Normal(int index) const;
    inline int GlobalIndex() const;
    inline bool IsNeighbors(const CurveEdge* other) const;
};

struct CurveVertex {
    PolyCurveNetwork* curve = nullptr;
    int id = 0;
    int component = 0;

    bool operator ==(const CurveVertex& other) const;
    bool operator !=(const CurveVertex& other) const;
    CurveVertex& operator +=(const Vector3& other);

    inline int numEdges() const;
    inline CurveEdge* edge(int i) const;
    inline CurveEdge* edge(size_t i) const;
    inline Vector3 Position() const;
    inline void SetPosition(const Vector3& pos) const;
    inline Vector3 Tangent() const;
    inline double DualLength() const;
    inline int GlobalIndex() const;
    inline Vector3 TotalLengthGradient() const;
    inline bool IsNeighbors(CurveVertex* other) const;
};

class PolyCurveNetwork {
    friend struct CurveVertex;

public:
    PolyCurveNetwork(std::vector<Vector3>& ps, std::vector<std::array<size_t, 2>>& es);
    PolyCurveNetwork(Eigen::MatrixXd& ps, std::vector<std::array<size_t, 2>>& es);
    ~PolyCurveNetwork();

    void InitStructs(std::vector<std::array<size_t, 2>>& es);
    void FindComponents();
    void PinVertex(int i);
    void PinTangent(int i);
    void PinToSurface(int i);
    void PinAllSpecialVertices(bool includeTangents);
    void PinAllEndpoints(bool includeTangents);
    void PrintPins() const;

    bool pinnedAllToSurface = false;

    inline Vector3 Position(const CurveVertex* v) const {
        return SelectRow(positions, v->id);
    }

    inline void SetPosition(const CurveVertex* v, const Vector3& newPos) {
        SetPosition(v->id, newPos);
    }

    std::vector<CurveVertex*> Vertices() const {
        return vertices;
    }

    std::vector<CurveEdge*> Edges() const {
        return edges;
    }

    inline std::vector<Vector3> AllPositions() const {
        std::vector<Vector3> nodes;
        for (const CurveVertex* vertex : vertices) {
            nodes.push_back(Position(vertex));
        }
        return nodes;
    }

    inline int NumVertices() const {
        return nVerts;
    }

    inline int NumEdges() const {
        return edges.size();
    }

    inline int NumPins() const {
        return pinnedVertices.size();
    }

    inline int NumTangentPins() const {
        return pinnedTangents.size();
    }

    inline int NumPinnedToSurface() const {
        return pinnedToSurface.size();
    }

    inline int NumComponents() const {
        return verticesByComponent.size();
    }

    inline int NumVerticesInComponent(int i) const {
        return verticesByComponent[i].size();
    }

    inline CurveVertex* GetVertexInComponent(int c, int i) const {
        return verticesByComponent[c][i];
    }

    inline CurveVertex* GetVertex(int i) const {
        return vertices[i];
    }

    inline CurveVertex* GetVertex(size_t i) const {
        return vertices[i];
    }

    inline CurveVertex* GetPinnedVertex(int i) const {
        return vertices[pinnedVertices[i]];
    }

    inline CurveVertex* GetPinnedTangent(int i) const {
        return vertices[pinnedTangents[i]];
    }

    inline CurveVertex* GetPinnedToSurface(int i) const {
        return vertices[pinnedToSurface[i]];
    }

    inline bool isPinned(int i) const {
        return pinnedSet.count(i) > 0;
    }

    inline bool isTangentPinned(int i) const {
        return pinnedTangentSet.count(i) > 0;
    }

    inline CurveEdge* GetEdge(int i) const {
        return edges[i];
    }

    inline CurveEdge* GetEdge(size_t i) const {
        return edges[i];
    }

    inline double BendingAngle(int i) const {
        CurveVertex* v_i = GetVertex(i);
        if (v_i->numEdges() != 2) return 0;

        CurveVertex* v_prev = v_i->edge(0)->Opposite(v_i);
        CurveVertex* v_next = v_i->edge(1)->Opposite(v_i);

        Vector3 e_prev = (v_i->Position() - v_prev->Position()).normalize();
        Vector3 e_next = (v_next->Position() - v_i->Position()).normalize();

        double sin_angle = norm(cross(e_prev, e_next));
        double cos_angle = dot(e_prev, e_next);
        return atan2(sin_angle, cos_angle);
    }

    void BoundingCube(Vector3& center, double& width) const;
    Vector3 Barycenter() const;
    double TotalLength() const;
    Vector3 AreaVector() const;
    PolyCurveNetwork* Subdivide() const;
    PolyCurveNetwork* Coarsen(MatrixProjectorOperator* op, bool doEdgeMatrix = false) const;

    NullSpaceProjector* constraintProjector;

    template <typename T>
    void AddConstraintProjector(DomainConstraints<T>& constraints) {
        if (constraintProjector) {
            delete constraintProjector;
        }
        constraintProjector = new NullSpaceProjector(constraints);
    }

    Eigen::MatrixXd positions;
    std::vector<ConstraintType> appliedConstraints;
    ImplicitSurface* constraintSurface;

private:
    int nVerts = 0;
    std::vector<int> pinnedToSurface;
    std::vector<int> pinnedVertices;
    std::unordered_set<int> pinnedSet;
    std::vector<int> pinnedTangents;
    std::unordered_set<int> pinnedTangentSet;
    std::vector<CurveVertex*> vertices;
    std::vector<CurveEdge*> edges;
    std::vector<std::vector<CurveEdge*>> adjacency;

    std::vector<std::vector<CurveVertex*>> verticesByComponent;
    std::vector<std::vector<CurveEdge*>> edgesByComponent;

    void CleanUpStructs();

    inline Vector3 Position(int i) const {
        return SelectRow(positions, i);
    }

    inline void SetPosition(int i, const Vector3& newPos) {
        SetRow(positions, i, newPos);
    }
};

inline CurveVertex* CurveEdge::Opposite(const CurveVertex* v) const {
    if (nextVert == v) return prevVert;
    else if (prevVert == v) return nextVert;
    else return 0;
}

inline Vector3 CurveEdge::Vector() const {
    return (nextVert->Position() - prevVert->Position());
}

inline Vector3 CurveEdge::Tangent() const {
    return this->Vector().normalize();
}

inline double CurveEdge::Length() const {
    return norm(nextVert->Position() - prevVert->Position());
}

inline Vector3 CurveEdge::Midpoint() const {
    return (nextVert->Position() + prevVert->Position()) / 2;
}

inline int CurveEdge::GlobalIndex() const {
    return id;
}

inline bool CurveEdge::IsNeighbors(const CurveEdge* other) const {
    const bool shared = (nextVert == other->nextVert) ||
                        (nextVert == other->prevVert) ||
                        (prevVert == other->nextVert) ||
                        (prevVert == other->prevVert);
    return shared && (other != this);
}

inline int CurveVertex::numEdges() const {
    return curve->adjacency[id].size();
}

inline CurveEdge* CurveVertex::edge(int i) const {
    return curve->adjacency[id][i];
}

inline CurveEdge* CurveVertex::edge(size_t i) const {
    return curve->adjacency[id][i];
}

inline Vector3 CurveVertex::Position() const {
    return SelectRow(curve->positions, id);
}

inline void CurveVertex::SetPosition(const Vector3& pos) const {
    curve->SetPosition(id, pos);
}

inline Vector3 CurveVertex::Tangent() const {
    Vector3 tangent{0, 0, 0};
    for (int i = 0; i < numEdges(); i++) {
        tangent += edge(i)->Tangent();
    }
    tangent = tangent.normalize();
    return tangent;
}

inline double CurveVertex::DualLength() const {
    double length = 0;
    for (int i = 0; i < numEdges(); i++) {
        length += edge(i)->Length();
    }
    return length / numEdges();
}

inline int CurveVertex::GlobalIndex() const {
    return id;
}

inline Vector3 CurveVertex::TotalLengthGradient() const {
    Vector3 total{0, 0, 0};
    // Moving vertex i only affects the lengths of the two edges
    // on either side of i.
    for (int e = 0; e < numEdges(); e++) {
        CurveEdge* e_i = edge(e);
        // Directions are away from the neighboring vertices.
        Vector3 inward = Position() - e_i->Opposite(this)->Position();
        // Magnitudes are 1, since length changes at 1 unit / unit.
        total += inward.normalize();
    }
    return total;
}

inline bool CurveVertex::IsNeighbors(CurveVertex* other) const {
    for (int i = 0; i < numEdges(); i++) {
        CurveEdge* e_i = edge(i);
        if (e_i->Opposite(this) == other) return true;
    }
    return false;
}
} // namespace LWS

