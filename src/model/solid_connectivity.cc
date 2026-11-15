#include "model/solid_connectivity.h"

#include <algorithm>
#include <vector>

#include "elements/solid_element_read_view.h"
namespace fesa {
namespace {
// Local primitive incidence record; no borrowed object survives this call.
struct Face {
  int element_id = 0;
  int corner_count = 0;
  int node_count = 0;
  int corners[4] = {};
  int nodes[6] = {};
};
bool Less(const Face* a, const Face* b) {
  if (a->corner_count != b->corner_count)
    return a->corner_count < b->corner_count;
  for (int i = 0; i < a->corner_count; ++i)
    if (a->corners[i] != b->corners[i]) return a->corners[i] < b->corners[i];
  return a->element_id < b->element_id;
}
bool SameCorners(const Face* a, const Face* b) {
  return a->corner_count == b->corner_count &&
         std::equal(a->corners, a->corners + a->corner_count, b->corners);
}
void Add(const SolidElementReadView* element, const int* indices, int corners,
         int nodes, std::vector<Face>* faces) {
  Face face;
  face.element_id = element->element()->Id();
  face.corner_count = corners;
  face.node_count = nodes;
  for (int i = 0; i < nodes; ++i)
    face.nodes[i] = element->GetNode(indices[i])->Id();
  for (int i = 0; i < corners; ++i) face.corners[i] = face.nodes[i];
  std::sort(face.corners, face.corners + corners);
  std::sort(face.nodes, face.nodes + nodes);
  faces->push_back(face);
}
StatusCode Unsupported(const Face* face, const char* reason, Diagnostic* d) {
  const auto status = ReportFailure(StatusCode::kUnsupported,
                                    DiagnosticStage::kValidation, reason, d);
  if (d) d->element_id = face->element_id;
  return status;
}
}  // namespace
StatusCode ValidateSolidConnectivity(const Domain* domain, Diagnostic* d) {
  bool has_solids = false;
  for (const auto& entry : *domain->Elements()) {
    SolidElementReadView view;
    if (SolidElementReadView::Bind(entry.second, &view, nullptr) ==
        StatusCode::kOk) {
      has_solids = true;
      break;
    }
  }
  if (!has_solids) return StatusCode::kOk;
  // MSVC Debug allocates an iterator proxy even for an empty vector. The
  // count constructor permits bad_alloc to reach the caller's boundary;
  // the default constructor is noexcept on this supported toolchain.
  std::vector<Face> faces(0);
  for (const auto& entry : *domain->Elements()) {
    SolidElementReadView view;
    const auto status = SolidElementReadView::Bind(entry.second, &view, nullptr);
    if (status == StatusCode::kUnsupported) continue;
    if (status != StatusCode::kOk) return status;
    const auto* element = &view;
    switch (element->formulation()) {
      case SolidFormulation::kTetra4Full:
      case SolidFormulation::kTetra10Full: {
        const int incidence[4][6] = {{0, 1, 2, 4, 5, 6},
                                     {0, 1, 3, 4, 8, 7},
                                     {1, 2, 3, 5, 9, 8},
                                     {0, 2, 3, 6, 9, 7}};
        for (const auto& face : incidence)
          Add(element, face, 3, element->node_count() == 10 ? 6 : 3, &faces);
        break;
      }
      case SolidFormulation::kHexa8Full:
      case SolidFormulation::kHexa8Reduced: {
        const int incidence[6][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
                                     {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
        for (const auto& face : incidence) Add(element, face, 4, 4, &faces);
        break;
      }
      case SolidFormulation::kWedge6: {
        const int triangles[2][3] = {{0, 1, 2}, {3, 4, 5}},
                  quads[3][4] = {{0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};
        for (const auto& face : triangles) Add(element, face, 3, 3, &faces);
        for (const auto& face : quads) Add(element, face, 4, 4, &faces);
        break;
      }
      case SolidFormulation::kPyramid5Full:
      case SolidFormulation::kPyramid5Reduced: {
        const int quad[] = {0, 1, 2, 3},
                  triangles[4][3] = {
                      {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};
        Add(element, quad, 4, 4, &faces);
        for (const auto& face : triangles) Add(element, face, 3, 3, &faces);
        break;
      }
    }
  }
  // Comparator references are the standard-library calling convention.
  const auto less = [](const Face& a, const Face& b) { return Less(&a, &b); };
  std::sort(faces.begin(), faces.end(), less);
  for (std::size_t i = 1; i < faces.size(); ++i) {
    const Face* a = &faces[i - 1];
    const Face* b = &faces[i];
    if (SameCorners(a, b) &&
        (a->node_count != b->node_count ||
         !std::equal(a->nodes, a->nodes + a->node_count, b->nodes)))
      return Unsupported(b, "Nonconforming solid face interpolation", d);
  }
  for (const auto& quad : faces) {
    if (quad.corner_count != 4) continue;
    for (int omit = 0; omit < 4; ++omit) {
      Face triangle;
      triangle.corner_count = 3;
      int count = 0;
      for (int corner = 0; corner < 4; ++corner)
        if (corner != omit) triangle.corners[count++] = quad.corners[corner];
      const auto found =
          std::lower_bound(faces.begin(), faces.end(), triangle, less);
      if (found != faces.end() && SameCorners(&*found, &triangle))
        return Unsupported(
            &quad, "Solid quad/triangle face requires conforming transition",
            d);
    }
  }
  return StatusCode::kOk;
}
}  // namespace fesa
