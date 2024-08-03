#include <future> 
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <iostream>
#include <limits>
#include <algorithm> 
#include <vector>
#include <string>
#include <cmath>
#include <regex>
#include <map>
#include <windows.h>
#include "rang.hpp"
#include "OCR_font_STL.h"

//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <boost/graph/graph_traits.hpp>
//#include <CGAL/boost/graph/IO/STL.h>

#include <vtkNew.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAutoInit.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkIdList.h>
#include <vtkAnnotatedCubeActor.h>
#include <vtkSphereSource.h>
#include <vtkAxesActor.h>
#include <vtkDiskSource.h>
#include <vtkCaptionActor2D.h>
#include <vtkSmartPointer.h>
#include <vtkTextProperty.h>
#include <vtkHardwarePicker.h>
#include <vtkTextActor.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkProperty2D.h>
#include <vtkSliderWidget.h>
#include <vtkCommand.h>
#include <vtkPlaneSource.h>
#include <vtkPlaneWidget.h>
#include <vtkTexturedButtonRepresentation2D.h>
#include <vtkButtonWidget.h>
#include <vtkImageData.h>
#include <vtkFreeTypeTools.h>
#include <vtkImageCanvasSource2D.h>
#include <vtkOpenGLRenderWindow.h>

#include <vtkWin32OpenGLRenderWindow.h>
#include <vtkSplineWidget.h>
#include <vtkSplineWidget2.h>
#include <vtkSplineRepresentation.h>
#include <vtkKochanekSpline.h>

#include <vtkButtonWidget.h>
#include <vtkCoordinate.h>
#include <vtkImageData.h>
#include <vtkNamedColors.h>
#include <vtkCallbackCommand.h>
#include <vtkCameraOrientationWidget.h>

#include <vtkDoubleArray.h>
#include <vtkCameraOrientationRepresentation.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCellData.h>
#include <vtkMath.h>

#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCellData.h>
#include <vtkMath.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkGlyph3D.h>
#include <vtkVertexGlyphFilter.h>



#define M_PI 3.14159265358979323846

//<< Red << << ColorEnd <<
auto ColorEnd = [](std::ostream& os) -> std::ostream& { return os << rang::fg::reset; };
auto Red = [](std::ostream& os) -> std::ostream& { return os << rang::fg::red; };
auto Green = [](std::ostream& os) -> std::ostream& { return os << rang::fg::green; };
auto Yellow = [](std::ostream& os) -> std::ostream& { return os << rang::fg::yellow; };
auto Blue = [](std::ostream& os) -> std::ostream& { return os << rang::fg::blue; };
auto Magenta = [](std::ostream& os) -> std::ostream& { return os << rang::fg::magenta; };
auto Cyan = [](std::ostream& os) -> std::ostream& { return os << rang::fg::cyan; };
auto Gray = [](std::ostream& os) -> std::ostream& { return os << rang::fg::gray; };

bool DEBUG = true;
namespace fs = std::filesystem;
namespace PMP = CGAL::Polygon_mesh_processing;

//typedef CGAL::Simple_cartesian<double>							Cgal_Kernel;
//typedef CGAL::Exact_predicates_exact_constructions_kernel		Cgal_Kernel;

typedef CGAL::Exact_predicates_inexact_constructions_kernel		Cgal_Kernel;
typedef CGAL::Surface_mesh<Cgal_Kernel::Point_3>				Cgal_Mesh;
//typedef Cgal_Kernel::FT											Cgal_ExactNumber;
typedef Cgal_Kernel::Point_3									Cgal_Point;
typedef Cgal_Kernel::Vector_3									Cgal_Vector;
typedef Cgal_Kernel::Plane_3									Cgal_Plane;
typedef Cgal_Kernel::Iso_cuboid_3								Cgal_cuboid;
typedef Cgal_Kernel::Aff_transformation_3						Cgal_Transformation;
typedef Cgal_Mesh::Vertex_index									Cgal_Vertex;
typedef Cgal_Mesh::Halfedge_index								Cgal_Halfedge;
typedef Cgal_Mesh::Face_index									Cgal_Face;
typedef boost::graph_traits<Cgal_Mesh>::face_descriptor			Cgal_face_descriptor;

VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
class C_InteractorStyle;

vtkNew<vtkActor> staticActor;
vtkNew<vtkActor> movableActor;


double maxcut = 10.0, mincut = -10.0;
int model_count = 0;

bool Duplicates = false;
int precision = 4;

struct PTSPoint {
	double x, y, z, u, v, w;
};

bool isXYZLine(const std::string& line) {
	std::istringstream iss(line);
	double x, y, z, u, v, w;
	iss >> x >> y >> z;
	if (iss.fail()) {
		return false; // Failed to parse XYZ
	}

	iss >> u >> v >> w; // Attempt to parse UVW
	if (!iss.fail() && iss.eof()) {
		return true; // Successfully parsed XYZUVW
	}

	iss.clear(); // Clear the fail state and check if we reached the end of the line after XYZ
	std::string remaining;
	iss >> remaining;
	return remaining.empty(); // True if only XYZ were in the line and nothing else
}

std::vector<PTSPoint> ReadPointsFromFile(const std::string& filename) {
	std::vector<PTSPoint> points;
	std::ifstream inFile(filename);
	std::string line;
	PTSPoint point;

	while (std::getline(inFile, line)) {
		if (!isXYZLine(line)) {
			continue; // Skip lines until an XYZ line is found
		}

		std::istringstream iss(line);
		if (iss >> point.x >> point.y >> point.z) { // Assuming the format is "x y z"
			points.push_back(point);
		}
	}

	return points;
}

bool needsPrecisionAdjustment(const std::string& coordinate) {
	size_t dotPosition = coordinate.find('.');
	if (dotPosition != std::string::npos && dotPosition + precision + 1 < coordinate.length()) {
		return true; // Check if there are more than x digits after the decimal point
	}
	return false;
}

std::string adjustPrecision(double value) {
	std::ostringstream oss;
	oss << std::fixed << std::setprecision(precision) << value;
	return oss.str();
}

std::string process_line(const std::string& line, double& offset_x, double& offset_y, double& offset_z
	, double angle_x, double angle_y, double angle_z) {
	std::istringstream iss(line);
	double x, y, z, temp_x, temp_y, temp_z,
		u = std::numeric_limits<double>::quiet_NaN(),
		v = std::numeric_limits<double>::quiet_NaN(),
		w = std::numeric_limits<double>::quiet_NaN();
	iss >> x >> y >> z;
	iss >> u >> v >> w; // Attempt to read UVW. Will fail for XYZ lines.
	bool uvwRead = !iss.fail(); // Check if UVW read was successful

	// Translate
	x -= offset_x;
	y -= offset_y;
	z -= offset_z;

	// Rotation around Z-axis
	temp_x = x * cos(angle_z * M_PI / 180.0) - y * sin(angle_z * M_PI / 180.0);
	temp_y = x * sin(angle_z * M_PI / 180.0) + y * cos(angle_z * M_PI / 180.0);
	x = temp_x;
	y = temp_y;

	// Rotation around Y-axis
	temp_x = x * cos(angle_y * M_PI / 180.0) + z * sin(angle_y * M_PI / 180.0);
	temp_z = -x * sin(angle_y * M_PI / 180.0) + z * cos(angle_y * M_PI / 180.0);
	x = temp_x;
	z = temp_z;

	// Rotation around X-axis
	temp_y = y * cos(angle_x * M_PI / 180.0) - z * sin(angle_x * M_PI / 180.0);
	temp_z = y * sin(angle_x * M_PI / 180.0) + z * cos(angle_x * M_PI / 180.0);
	y = temp_y;
	z = temp_z;

	std::string sx = std::to_string(x);
	std::string sy = std::to_string(y);
	std::string sz = std::to_string(z);
	std::string su, sv, sw;


	if (uvwRead) {
		su = std::to_string(u);
		sv = std::to_string(v);
		sw = std::to_string(w);
	}

	bool needsAdjustment = needsPrecisionAdjustment(sx) || needsPrecisionAdjustment(sy) || needsPrecisionAdjustment(sz) ||
		(uvwRead && (needsPrecisionAdjustment(su) || needsPrecisionAdjustment(sv) || needsPrecisionAdjustment(sw)));

	if (needsAdjustment) {
		std::ostringstream oss;
		oss << adjustPrecision(x) << " "
			<< adjustPrecision(y) << " "
			<< adjustPrecision(z);
		if (uvwRead) {
			oss << " " << adjustPrecision(u) << " " << adjustPrecision(v) << " " << adjustPrecision(w);
		}
		return oss.str();
	}
	else {
		return line;
	}
}

std::vector<std::string> removeConsecutiveDuplicates(const std::vector<std::string>& lines) {
	std::vector<std::string> uniqueLines;
	std::unordered_set<std::string> seen;

	for (const auto& line : lines) {
		// Attempt to insert the line into the set to check if it's a duplicate
		if (seen.insert(line).second) { // insert returns a pair, where .second is true if the insertion took place
			uniqueLines.push_back(line); // If the line was inserted successfully, it's not a duplicate
		}
		else {
			Duplicates = true;
		}
	}

	return uniqueLines;
}

void process_pts(std::string& inputFile, std::string& outputFile, double& offset_x, double& offset_y, double& offset_z,
	double& rotX, double& rotY, double& rotZ) {
	bool _ReadLine = true;
	std::vector<std::string> lines;
	std::ifstream PTS_in(inputFile);
	if (PTS_in) {
		std::vector<PTSPoint> points = ReadPointsFromFile(inputFile);
		std::string line;

		while (std::getline(PTS_in, line)) {
			if (!isXYZLine(line)) {
				continue; // Skip lines until an XYZ line is found
			}
			line = process_line(line, offset_x, offset_y, offset_z, rotX, rotY, rotZ); // Process the line to fixed precision
			lines.push_back(line);
		}
		lines = removeConsecutiveDuplicates(lines); // Remove duplicates while preserving order
	}

	std::ofstream PTS_out(outputFile);
	if (PTS_out) {
		PTS_out << "# (Created by Banna)" << "\n";
		for (const auto& l : lines) {
			PTS_out << l << "\n";
		}
	}
}

static void displayUserName() {
	char* username = nullptr;
	char* userdomain = nullptr;
	size_t sizeUsername = 0;
	size_t sizeUserdomain = 0;

	_dupenv_s(&username, &sizeUsername, "USERNAME");
	_dupenv_s(&userdomain, &sizeUserdomain, "USERDOMAIN");

	std::cout << "\n        USERNAME: "
		<< (userdomain ? userdomain : "Unknown") << "\\"
		<< (username ? username : "Unknown") << std::endl;

	free(username);
	free(userdomain);
}


void setConsoleSize(int width, int height) {
	HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); // Get the standard output handle

	COORD newSize;
	newSize.X = width;
	newSize.Y = 32766; // Maximal possible height for the console window
	SetConsoleScreenBufferSize(hStdout, newSize);

	SMALL_RECT windowSize;
	windowSize.Top = 0;
	windowSize.Left = 0;
	windowSize.Right = width - 1;  // Width of the window
	windowSize.Bottom = height - 1;  // Height of the window

	if (!SetConsoleWindowInfo(hStdout, TRUE, &windowSize)) {
		std::cerr << "        Setting console window size failed." << std::endl;
	}
}

class RotationSliderCallback : public vtkCommand {
public:
	static RotationSliderCallback* New() {
		return new RotationSliderCallback;
	}

	RotationSliderCallback() : RotActor(nullptr) {}

	void Execute(vtkObject* caller, unsigned long, void*) override {
		vtkSliderWidget* sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
		vtkSliderRepresentation* sliderRep = dynamic_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation());
		if (!sliderRep) return;
		double value = scaleValue(sliderRep->GetValue(), 180.0);

		if (RotPtr) *RotPtr = value;  // Update the external CutHeight variable
		else std::cerr << "        Warning: RotPtr is not initialized." << std::endl;

		if (this->RotActor) {
			double CurrentRot[4];
			this->RotActor->GetOrientation(CurrentRot);
			this->RotActor->SetOrientation(CurrentRot[0], CurrentRot[1], value);  // Set Z-position of the cutting plane
			//if (DEBUG) std::cout << Yellow << "        Rotation Slider : " << ColorEnd << value << std::endl;
			//if (DEBUG) std::cout << Yellow << "        Mesh Orientation:  " << ColorEnd
			//	<< CurrentRot[0] << "  "
			//	<< CurrentRot[1] << "  "
			//	<< CurrentRot[2] << "  "
			//	<< CurrentRot[3] << std::endl;
		}
		char label[50];
		sprintf_s(label, "%.1f", value);  // Format to two decimal places
		sliderRep->SetLabelFormat(label);

		sliderWidget->GetInteractor()->GetRenderWindow()->Render(); // Update the display
	}

	void SetRotActor(vtkActor* actor) {
		this->RotActor = actor;
	}
	void SetRotPtr(double* ptr) {
		RotPtr = ptr;
	}
	double scaleValue(double input, double factor) {
		double normalizedInput = input / factor;  // Normalize to -1 to 1
		return factor * normalizedInput * normalizedInput * (input < 0 ? -1 : 1);  // Scale back to -180 to 180
	}

private:
	vtkActor* RotActor;
	double* RotPtr = nullptr;
};

class CutSliderCallback : public vtkCommand
{
public:
	static CutSliderCallback* New() {
		return new CutSliderCallback();
	}

	CutSliderCallback() : CutHeightPtr(nullptr) {}

	virtual void Execute(vtkObject* caller, unsigned long eventId, void*) override {
		vtkSliderWidget* sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
		vtkSliderRepresentation* sliderRep = dynamic_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation());
		if (!sliderRep) return;
		double value = sliderRep->GetValue();
		if (CutHeightPtr) *CutHeightPtr = value;  // Update the external CutHeight variable
		else std::cerr << "        Warning: CutHeightPtr is not initialized." << std::endl;

		if (value >= mincut || value <= maxcut) {
			double CurrentPos[4];
			this->CuttingDisk->GetPosition(CurrentPos);
			this->CuttingDisk->SetPosition(CurrentPos[0], CurrentPos[1], value);  // Set Z-position of the cutting plane
			//if (DEBUG) std::cout << Yellow << "        Cut Slider : " << ColorEnd << value << std::endl;
		}

		char label[50];
		sprintf_s(label, "%.1f", value);  // Format to two decimal places
		sliderRep->SetLabelFormat(label);

		sliderWidget->GetInteractor()->GetRenderWindow()->Render(); // Update the display
	}

	void SetModelActor(vtkActor* actor) {
		this->CuttingDisk = actor;
	}
	void SetCutHeightPtr(double* ptr) {
		CutHeightPtr = ptr;
	}

private:
	vtkActor* CuttingDisk;
	double* CutHeightPtr = nullptr;
};

class C_InteractorStyle : public vtkInteractorStyleTrackballCamera {
public:
	static C_InteractorStyle* New();
	//vtkTypeMacro(C_InteractorStyle, vtkInteractorStyleTrackballCamera);

	vtkSmartPointer<vtkHardwarePicker> Picker;
	vtkSmartPointer<vtkTextActor> TextActor;
	vtkSmartPointer<vtkActor> SelectedMesh;  // The mesh to pick
	vtkSmartPointer<vtkActor> MeshActor;
	bool IsMeshSelected = false;
	int LastPosition[2] = { -1, -1 };
	double X_offset, Y_offset;

	C_InteractorStyle() : X_offset(0), Y_offset(0) {
		this->Picker = vtkSmartPointer<vtkHardwarePicker>::New();
		this->TextActor = vtkSmartPointer<vtkTextActor>::New();
		this->TextActor->GetTextProperty()->SetFontSize(20);
		this->TextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // 0.7, 0.5, 0.3
		this->TextActor->SetDisplayPosition(10, 10);
		this->TextActor->SetInput("Adjust the Model");

		this->SelectedMesh = nullptr;
		this->IsMeshSelected = false;
		this->LastPosition[0] = this->LastPosition[1] = 0;
	}

	virtual void OnLeftButtonDown() override
	{
		int* clickPos = this->GetInteractor()->GetEventPosition();
		this->Picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());
		vtkActor* Mactor = vtkActor::SafeDownCast(this->Picker->GetActor());

		if (Mactor && Mactor == this->MeshActor) {
			this->SelectedMesh = Mactor;
			this->IsMeshSelected = true;
			this->LastPosition[0] = clickPos[0];
			this->LastPosition[1] = clickPos[1];

			if (DEBUG) std::cout << Yellow << "        Mesh selected!" << ColorEnd << std::endl;
			this->TextActor->SetInput("Mesh selected!");
			this->SelectedMesh->GetProperty()->SetColor(1.0, 1.0, 0.0);  // Change color when selected
		}
		else {
			if (DEBUG) std::cout << Yellow << "        Nothing selected!" << std::endl;
			this->LastPosition[0] = 0; // Reset last position if no valid selection
			this->LastPosition[1] = 0;
		}
		this->Interactor->GetRenderWindow()->Render();
	}

	virtual void OnLeftButtonUp() override {
		if (this->IsMeshSelected && this->SelectedMesh) {
			double CurrentPos[4];
			this->SelectedMesh->GetPosition(CurrentPos);

			X_offset = CurrentPos[0];
			Y_offset = CurrentPos[1];

			if (DEBUG) std::cout << Yellow << "        Final position offset: "
				<< ColorEnd << "X " << CurrentPos[0] << "  Y " << CurrentPos[1] << std::endl;

			this->TextActor->SetInput(("X " + std::format("{:.2f}", CurrentPos[0]) + " Y " + std::format("{:.2f}", CurrentPos[1])).c_str());
			this->SelectedMesh->GetProperty()->SetColor(0.85, 0.85, 0.85);
			this->Interactor->GetRenderWindow()->Render();

			this->IsMeshSelected = false;
			this->SelectedMesh = nullptr;
		}
	}

	virtual void OnMouseMove() override {
		if (this->IsMeshSelected && this->SelectedMesh) {
			vtkRenderWindowInteractor* RenderWindowInteractor = this->Interactor;
			vtkRenderer* renderer = this->GetDefaultRenderer();

			int* newPos = RenderWindowInteractor->GetEventPosition();

			// Convert new mouse position to world coordinates
			double displayPos[3] = { double(newPos[0]), double(newPos[1]), 0.0 };
			renderer->SetDisplayPoint(displayPos);
			renderer->DisplayToWorld();
			renderer->GetWorldPoint(displayPos);

			// Convert last position to world coordinates
			double lastDisplayPos[3] = { double(LastPosition[0]), double(LastPosition[1]), 0.0 };
			renderer->SetDisplayPoint(lastDisplayPos);
			renderer->DisplayToWorld();
			renderer->GetWorldPoint(lastDisplayPos);

			// Calculate movement deltas
			double Dx = displayPos[0] - lastDisplayPos[0];
			double Dy = displayPos[1] - lastDisplayPos[1];

			// Update mesh position
			double pos[3];
			this->SelectedMesh->GetPosition(pos);
			this->SelectedMesh->SetPosition(pos[0] + Dx, pos[1] + Dy, pos[2]);

			//if (DEBUG) std::cout << Yellow << "        Moved to " << ColorEnd
			//	<< "X " << std::format("{:.2f}", pos[0] + Dx) << ", Y " << std::format("{:.2f}", pos[1] + Dy) << std::endl;
			this->TextActor->SetInput(("X " + std::format("{:.2f}", pos[0] + Dx) + " Y " + std::format("{:.2f}", pos[1] + Dy)).c_str());

			this->LastPosition[0] = newPos[0];
			this->LastPosition[1] = newPos[1];

			RenderWindowInteractor->Render();
		}

		if (!this->IsMeshSelected) {
			vtkInteractorStyleTrackballCamera::OnMouseMove();
		}
	}

	virtual void OnMiddleButtonDown() override {
		this->StartPan();
	}

	virtual void OnMiddleButtonUp() override {
		this->EndPan();
	}

	virtual void OnRightButtonDown() override {
		this->StartRotate();
	}

	virtual void OnRightButtonUp() override {
		this->EndRotate();
	}

	virtual void Pan() override {
		if (this->CurrentRenderer == nullptr || this->Interactor == nullptr) return;
		vtkRenderWindowInteractor* rwi = this->Interactor;
		vtkCamera* camera = this->CurrentRenderer->GetActiveCamera();
		if (!camera) return;
		int* lastPos = rwi->GetLastEventPosition();
		int* newPos = rwi->GetEventPosition();
		double dx = newPos[0] - lastPos[0];
		double dy = newPos[1] - lastPos[1];
		double scale = 0.05; // Adjust this scale to control the sensitivity of panning
		dx *= scale;
		dy *= scale;
		double right[3], up[3];
		camera->GetViewUp(up);
		this->CurrentRenderer->GetActiveCamera()->OrthogonalizeViewUp();
		vtkMath::Cross(camera->GetDirectionOfProjection(), up, right);
		vtkMath::Normalize(right);
		double cameraPosition[3], cameraFocalPoint[3];
		camera->GetPosition(cameraPosition);
		camera->GetFocalPoint(cameraFocalPoint);
		for (int i = 0; i < 3; i++) {
			cameraPosition[i] += dx * right[i] + dy * up[i];
			cameraFocalPoint[i] += dx * right[i] + dy * up[i];
		}
		camera->SetPosition(cameraPosition);
		camera->SetFocalPoint(cameraFocalPoint);
		this->CurrentRenderer->ResetCameraClippingRange();
		rwi->Render();
	}

	virtual void Dolly(double amount) override {
		if (this->CurrentRenderer == nullptr || this->Interactor == nullptr) return;
		vtkCamera* camera = this->CurrentRenderer->GetActiveCamera();
		if (!camera) return;
		const double DollyScaleFactor = 0.1, minDis = 50, maxDis = 300;
		amount = 1.0 + (amount - 1.0) * DollyScaleFactor;

		double newDistance = (camera->GetDistance()) * (1.0 / amount);  // Adjust the interpretation of amount

		if (newDistance < minDis) camera->SetDistance(minDis);
		else if (newDistance > maxDis) camera->SetDistance(maxDis);
		else camera->Dolly(amount);

		//if (DEBUG) std::cout << Yellow << "        Current Zoom Level: " << ColorEnd << camera->GetDistance() << std::endl;
		this->CurrentRenderer->ResetCameraClippingRange();
		this->Interactor->Render();
	}
};
vtkStandardNewMacro(C_InteractorStyle);

vtkNew<vtkPolyData> mesh_to_vtk(const Cgal_Mesh& mesh, bool DEBUG) {
	vtkNew<vtkPoints> VTKpoints;
	vtkNew<vtkCellArray> polygons;
	std::map<Cgal_Vertex, vtkIdType> vertexIdMap;
	vtkIdType vtkId = 0;
	for (auto v : mesh.vertices()) {
		const Cgal_Point& p = mesh.point(v);
		//VTKpoints->InsertNextPoint(CGAL::to_double(p.x()), CGAL::to_double(p.y()),  CGAL::to_double(p.z()));
		VTKpoints->InsertNextPoint(p.x(), p.y(), p.z());
		vertexIdMap[v] = vtkId++;
	}
	for (auto f : mesh.faces()) {
		vtkNew<vtkIdList> polygon;
		CGAL::Vertex_around_face_iterator<Cgal_Mesh> vbegin, vend;
		boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(f), mesh);
		for (; vbegin != vend; ++vbegin) {
			Cgal_Vertex v = *vbegin;
			polygon->InsertNextId(vertexIdMap[v]);
			//polygon->InsertNextId(*vbegin);
		}
		polygons->InsertNextCell(polygon);
	}
	vtkNew<vtkPolyData> polyData;
	polyData->SetPoints(VTKpoints);
	polyData->SetPolys(polygons);
	if (DEBUG) std::cout << Yellow << "        Mesh prepared for Viewer." << ColorEnd << std::endl;
	return polyData;
}

void CreateBTNImage(vtkImageData* image, int Xn, int Yn, unsigned char* color) {
	image->SetDimensions(Xn, Yn, 1);
	image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	for (int y = 0; y < Yn; y++) {
		for (int x = 0; x < Xn; x++) {
			unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
			for (int i = 0; i < 3; ++i) {
				pixel[i] = color[i];
			}
		}
	}
}

void CenterTextActor(vtkTextActor* textActor, vtkRenderer* renderer, double buttonMidX, double buttonMidY) {
	size_t CharN = strlen(textActor->GetInput());
	double textPosX = buttonMidX - (textActor->GetTextProperty()->GetFontSize() * 0.375 * CharN),
		textPosY = buttonMidY - (textActor->GetTextProperty()->GetFontSize() * 0.6);
	textActor->SetPosition(textPosX, textPosY);
}


void OKButtonCallbackFunction(vtkObject* caller, long unsigned int, void* clientData, void*) {
	vtkTextActor* textActor = static_cast<vtkTextActor*>(clientData);
	if (textActor == nullptr) return;
	const char* currentText = textActor->GetInput();
	if (strcmp(currentText, "OK") == 0) {
		textActor->SetInput("DONE");
		std::cout << "DONE" << std::endl;
	}
	else {
		textActor->SetInput("OK");
		std::cout << "OK" << std::endl;
	}

	// Get the button representation and renderer from the button widget
	vtkButtonWidget* buttonWidget = static_cast<vtkButtonWidget*>(caller);
	vtkTexturedButtonRepresentation2D* buttonRep = static_cast<vtkTexturedButtonRepresentation2D*>(buttonWidget->GetRepresentation());
	vtkRenderer* renderer = buttonRep->GetRenderer();
	CenterTextActor(textActor, renderer, 60, 105);
	renderer->GetRenderWindow()->Render();
}

void CLOSEButtonCallbackFunction(vtkObject* caller, long unsigned int, void*, void*) {
	vtkButtonWidget* buttonWidget = static_cast<vtkButtonWidget*>(caller);
	vtkTexturedButtonRepresentation2D* buttonRep = static_cast<vtkTexturedButtonRepresentation2D*>(buttonWidget->GetRepresentation());
	vtkRenderer* renderer = buttonRep->GetRenderer();
	vtkRenderWindow* renderWindow = renderer->GetRenderWindow();
	vtkRenderWindowInteractor* interactor = renderWindow->GetInteractor();

	// Finalize and close the render window
	renderWindow->Finalize();
	interactor->TerminateApp();
	return;
}

//void MouseMoveCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData) {
//	vtkRenderWindowInteractor* interactor = static_cast<vtkRenderWindowInteractor*>(caller);
//	int x, y;
//	interactor->GetEventPosition(x, y);
//	std::cout << "        Mouse position: (" << x << ", " << y << ")" << std::endl;
//}

void visualize_mesh(Cgal_Mesh movableMesh, double Height_, const std::string ID, Cgal_Mesh staticMesh, double& Xoffset, double& Yoffset,
	double& CutHeight, double& RotZ, bool DEBUG) {

	if (DEBUG) std::cout << Yellow << "        Starting Mesh Viewer." << ColorEnd << std::endl;
	CutHeight = 0.0; RotZ = 0.0;

	byte WIN_TRANS = 245;
	double bkR = 0.129, bkG = 0.129, bkB = 0.141,
		movR = 1.0, movG = 1.0, movB = 1.0,
		stR = 0.7, stG = 0.5, stB = 0.3;
	//double bkR = 0.2, bkG = 0.2, bkB = 0.3;
	// Main Renderer setup
	vtkNew<vtkRenderer> MainRenderer;
	MainRenderer->SetBackground(bkR, bkG, bkB);
	MainRenderer->SetUseDepthPeeling(1);
	MainRenderer->SetMaximumNumberOfPeels(100);
	MainRenderer->SetOcclusionRatio(0.1);

	// Setup for secondary renderer
	vtkNew<vtkRenderer> InsetRenderer;
	InsetRenderer->SetBackground(bkR, bkG, bkB);
	InsetRenderer->SetUseDepthPeeling(1);
	InsetRenderer->SetMaximumNumberOfPeels(100);
	InsetRenderer->SetOcclusionRatio(0.1);
	InsetRenderer->SetViewport(0.75, 0.0, 1.0, 0.25); // X1 Y1   X2 Y2
	InsetRenderer->EraseOff();

	// Window and interactor setup
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->SetSize(700, 600);
	renderWindow->SetWindowName("AB Dental Tool");
	renderWindow->SetAlphaBitPlanes(1);
	renderWindow->SetMultiSamples(0);
	renderWindow->AddRenderer(MainRenderer);
	renderWindow->AddRenderer(InsetRenderer);

	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Mesh mappers and actors
	vtkNew<vtkPolyDataMapper> staticMapper;
	staticMapper->SetInputData(mesh_to_vtk(staticMesh, DEBUG));
	vtkNew<vtkActor> staticActor;
	staticActor->SetMapper(staticMapper);
	//staticActor->GetProperty()->SetColor(stR, stG, stB);
	staticActor->GetProperty()->SetDiffuse(0.9);
	staticActor->GetProperty()->SetDiffuseColor(stR, stG, stB);
	staticActor->GetProperty()->SetSpecular(0.3);
	staticActor->GetProperty()->SetSpecularPower(70.0);
	staticActor->GetProperty()->SetInterpolationToPhong();
	MainRenderer->AddActor(staticActor);

	vtkNew<vtkPolyDataMapper> movableMapper;
	movableMapper->SetInputData(mesh_to_vtk(movableMesh, DEBUG));
	vtkNew<vtkActor> movableActor;
	movableActor->SetMapper(movableMapper);
	//movableActor->GetProperty()->SetColor(movR, movG, movB);
	movableActor->GetProperty()->SetDiffuse(0.9);
	movableActor->GetProperty()->SetDiffuseColor(movR, movG, movB);
	movableActor->GetProperty()->SetSpecular(0.3);
	movableActor->GetProperty()->SetSpecularPower(70.0);
	movableActor->GetProperty()->SetInterpolationToPhong();
	MainRenderer->AddActor(movableActor);
	InsetRenderer->AddActor(movableActor);

	//// PTS actor
	//std::ifstream ptsFile(PTSfilepath);
	//vtkNew<vtkPoints> PTSpoints;
	//std::string line;
	//while (std::getline(ptsFile, line)) {
	//	std::istringstream iss(line);
	//	double x, y, z;
	//	iss >> x >> y >> z;
	//	PTSpoints->InsertNextPoint(x, y, z);
	//}
	//
	//vtkNew<vtkPolyData> PTSpointPolyData;
	//PTSpointPolyData->SetPoints(PTSpoints);
	//vtkNew<vtkVertexGlyphFilter> vertexFilter;
	//vertexFilter->SetInputData(PTSpointPolyData);
	//vertexFilter->Update();
	//vtkNew<vtkSphereSource> sphereSource;
	//sphereSource->SetRadius(0.75);  // Set the radius of the spheres
	//vtkNew<vtkGlyph3D> glyph3D;
	//glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
	//glyph3D->SetInputConnection(vertexFilter->GetOutputPort());
	//glyph3D->ScalingOff();  // Turn off scaling to keep spheres uniform
	//vtkNew<vtkPolyDataMapper> glyphMapper;
	//glyphMapper->SetInputConnection(glyph3D->GetOutputPort());
	//vtkNew<vtkActor> glyphActor;
	//glyphActor->SetMapper(glyphMapper);
	//glyphActor->GetProperty()->SetColor(stR, stG, stB);
	//MainRenderer->AddActor(glyphActor);

	// Model Name Text
	vtkNew<vtkTextActor> ModelTextActor;
	ModelTextActor->GetTextProperty()->SetFontSize(20);
	ModelTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // 0.7, 0.5, 0.3
	ModelTextActor->SetDisplayPosition(10, 35);
	ModelTextActor->SetInput(("ID: " + ID).c_str());
	MainRenderer->AddActor(ModelTextActor);

	// Model Height Text
	vtkNew<vtkTextActor> HeightTextActor;
	HeightTextActor->GetTextProperty()->SetFontSize(20);
	HeightTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // 0.7, 0.5, 0.3
	HeightTextActor->SetDisplayPosition(10, 60);
	HeightTextActor->SetInput(("Height: " + std::format("{:.2f}", Height_) + " mm").c_str());
	MainRenderer->AddActor(HeightTextActor);

	// Orientation Marker Widget - SetupCubeWidget(renderWindowInteractor, MainRenderer);
	vtkNew<vtkAnnotatedCubeActor> cubeActor;
	cubeActor->SetFaceTextScale(0.5);
	cubeActor->SetXPlusFaceText("X+");
	cubeActor->SetXMinusFaceText("X-");
	cubeActor->SetYPlusFaceText("Y+");
	cubeActor->SetYMinusFaceText("Y-");
	cubeActor->SetZPlusFaceText("+Z");
	cubeActor->SetZMinusFaceText("Z-");
	cubeActor->GetXPlusFaceProperty()->SetColor(0.7, 0.5, 0.3);   // Red for X+
	cubeActor->GetXMinusFaceProperty()->SetColor(0.7, 0.5, 0.3);  // Red for X-
	cubeActor->GetYPlusFaceProperty()->SetColor(0.7, 0.5, 0.3);   // Green for Y+
	cubeActor->GetYMinusFaceProperty()->SetColor(0.7, 0.5, 0.3);  // Green for Y-
	cubeActor->GetZPlusFaceProperty()->SetColor(0.7, 0.5, 0.3);   // Blue for Z+
	cubeActor->GetZMinusFaceProperty()->SetColor(0.7, 0.5, 0.3);  // Blue for Z-
	cubeActor->SetXFaceTextRotation(0);
	cubeActor->SetYFaceTextRotation(0);
	cubeActor->SetZFaceTextRotation(-90);
	cubeActor->GetCubeProperty()->SetColor(movR, movG, movB);
	vtkNew<vtkOrientationMarkerWidget> CubeActorWidget;
	CubeActorWidget->SetOrientationMarker(cubeActor);
	CubeActorWidget->SetViewport(0.0, 0.85, 0.15, 1.0);
	CubeActorWidget->SetInteractor(renderWindowInteractor);
	CubeActorWidget->SetDefaultRenderer(MainRenderer);
	CubeActorWidget->SetEnabled(1);
	CubeActorWidget->InteractiveOff();

	// Center XYZ axes - SetupCenterXYZAxes(renderWindowInteractor, MainRenderer);
	vtkNew<vtkAxesActor> XYZaxes;
	XYZaxes->SetTotalLength(6, 6, 6);
	XYZaxes->SetXAxisLabelText("");
	XYZaxes->SetYAxisLabelText("");
	XYZaxes->SetZAxisLabelText("");
	MainRenderer->AddActor(XYZaxes);

	// Center sphere - SetupSphereActor(MainRenderer);
	vtkNew<vtkSphereSource> sphereSource;
	sphereSource->SetCenter(0.0, 0.0, 0.0);
	sphereSource->SetRadius(0.5);
	sphereSource->Update();
	vtkNew<vtkPolyDataMapper> sphereMapper;
	sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
	vtkNew<vtkActor> sphereActor;
	sphereActor->SetMapper(sphereMapper);
	sphereActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
	MainRenderer->AddActor(sphereActor);

	// Create a center plane
	vtkNew<vtkPlaneSource> planeSource;
	planeSource->SetXResolution(10);
	planeSource->SetYResolution(10);
	planeSource->SetOrigin(-50, -50, 0);
	planeSource->SetPoint1(50, -50, 0);
	planeSource->SetPoint2(-50, 50, 0);
	vtkNew<vtkPolyDataMapper> planeMapper;
	planeMapper->SetInputConnection(planeSource->GetOutputPort());
	vtkNew<vtkActor> planeActor;
	planeActor->SetMapper(planeMapper);
	planeActor->GetProperty()->SetRepresentationToWireframe();
	planeActor->GetProperty()->SetColor(0.3, 0.3, 0.3);
	MainRenderer->AddActor(planeActor);

	// XY Cutting disk - SetupCutdiskActor(InsetRenderer);
	vtkNew<vtkDiskSource> diskSource;
	diskSource->SetInnerRadius(0.0);
	diskSource->SetOuterRadius(40.0);
	diskSource->SetRadialResolution(50);
	diskSource->SetCircumferentialResolution(50);
	diskSource->SetNormal(0, 0, 1);
	diskSource->Update();
	vtkNew<vtkPolyDataMapper> mapper;
	mapper->SetInputConnection(diskSource->GetOutputPort());
	vtkNew<vtkActor> diskActor;
	diskActor->SetMapper(mapper);
	//diskActor->GetProperty()->SetColor(0.2, 0.5, 0.4);
	diskActor->GetProperty()->SetColor(0.2, 0.2, 0.2);
	diskActor->GetProperty()->SetOpacity(0.5);
	InsetRenderer->AddActor(diskActor);

	// Cutting Slider - SetupCutSliderWidget(renderWindowInteractor, movableActor, CutHeight, mincut, maxcut, mincut);
	vtkNew<vtkSliderRepresentation2D> CutSlider;
	CutSlider->SetMinimumValue(mincut);
	CutSlider->SetMaximumValue(maxcut);
	CutSlider->SetValue(0.0);
	CutSlider->GetSliderProperty()->SetColor(0.7, 0.5, 0.3);
	CutSlider->GetTitleProperty()->SetColor(1.0, 1.0, 1.0);
	CutSlider->SetTubeWidth(0.005);
	CutSlider->SetSliderLength(0.05);
	CutSlider->SetSliderWidth(0.02);
	CutSlider->SetEndCapLength(0.005);
	CutSlider->SetEndCapWidth(0.02);
	CutSlider->SetTitleText("Cut");
	CutSlider->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
	CutSlider->GetPoint1Coordinate()->SetValue(0.9, 0.20);  // Bottom-left 
	CutSlider->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
	CutSlider->GetPoint2Coordinate()->SetValue(0.9, 0.80);  // Top-left
	vtkNew<vtkSliderWidget> CutSliderWidget;
	CutSliderWidget->SetInteractor(renderWindowInteractor);
	CutSliderWidget->SetRepresentation(CutSlider);
	CutSliderWidget->SetAnimationModeToAnimate();
	CutSliderWidget->SetEnabled(1);
	vtkNew<CutSliderCallback> CutSlidercallback;
	CutSlidercallback->SetModelActor(movableActor);
	CutSlidercallback->SetCutHeightPtr(&CutHeight);
	CutSliderWidget->AddObserver(vtkCommand::InteractionEvent, CutSlidercallback);

	// Rotate slider - SetupRotSliderWidget(renderWindowInteractor, movableActor, RotZ, -180, 180, 0);
	vtkNew<vtkSliderRepresentation2D> RotSlider;
	RotSlider->SetMinimumValue(-180);
	RotSlider->SetMaximumValue(180);
	RotSlider->SetValue(0);
	RotSlider->GetSliderProperty()->SetColor(0.7, 0.5, 0.3);
	RotSlider->GetTitleProperty()->SetColor(1.0, 1.0, 1.0);
	RotSlider->SetTubeWidth(0.005);
	RotSlider->SetSliderLength(0.05);
	RotSlider->SetSliderWidth(0.02);
	RotSlider->SetEndCapLength(0.005);
	RotSlider->SetEndCapWidth(0.02);
	RotSlider->SetTitleText("Rotate");
	RotSlider->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
	RotSlider->GetPoint1Coordinate()->SetValue(0.3, 0.1);
	RotSlider->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
	RotSlider->GetPoint2Coordinate()->SetValue(0.7, 0.1);
	vtkNew<vtkSliderWidget> RotSliderWidget;
	RotSliderWidget->SetInteractor(renderWindowInteractor);
	RotSliderWidget->SetRepresentation(RotSlider);
	RotSliderWidget->SetAnimationModeToAnimate();
	RotSliderWidget->SetEnabled(1);
	vtkNew<RotationSliderCallback> RotSlidercallback;
	RotSlidercallback->SetRotActor(movableActor);
	RotSlidercallback->SetRotPtr(&RotZ);
	RotSliderWidget->AddObserver(vtkCommand::InteractionEvent, RotSlidercallback);


	// Button Widget
	double width = 50, height = 15;
	unsigned char TextureColor[3] = { 179, 128, 77 };

	// OK Button Widget
	vtkNew<vtkImageData> Texture1;
	CreateBTNImage(Texture1, width, height, TextureColor);
	vtkNew<vtkTexturedButtonRepresentation2D> OKbuttonRepresentation;
	OKbuttonRepresentation->SetNumberOfStates(1);
	OKbuttonRepresentation->SetButtonTexture(0, Texture1);
	vtkNew<vtkButtonWidget> OKbuttonWidget;
	OKbuttonWidget->SetInteractor(renderWindowInteractor);
	OKbuttonWidget->SetRepresentation(OKbuttonRepresentation);
	double BTN1x = 60, BTN1y = 105;
	double bds1[6] = {
		BTN1x - width, BTN1x + width,  // Adjust width boundaries
		BTN1y - height, BTN1y + height,  // Adjust height boundaries
		0.0, 0.0
	};
	OKbuttonRepresentation->SetPlaceFactor(1);
	OKbuttonRepresentation->PlaceWidget(bds1);
	OKbuttonWidget->On();
	vtkNew<vtkTextActor> OKbtntextActor;
	OKbtntextActor->SetInput("OK");
	OKbtntextActor->GetTextProperty()->SetFontSize(20);
	OKbtntextActor->GetTextProperty()->SetBold(1);
	OKbtntextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
	CenterTextActor(OKbtntextActor, MainRenderer, BTN1x, BTN1y);
	MainRenderer->AddActor2D(OKbtntextActor);
	vtkNew<vtkCallbackCommand> OKbuttoncallback;
	OKbuttoncallback->SetCallback(OKButtonCallbackFunction);
	OKbuttoncallback->SetClientData(OKbtntextActor);
	OKbuttonWidget->AddObserver(vtkCommand::StateChangedEvent, OKbuttoncallback);

	// Exit Button Widget
	double BTN2x = 60, BTN2y = 155;
	vtkNew<vtkImageData> Texture2;
	CreateBTNImage(Texture2, width, height, TextureColor);
	vtkNew<vtkTexturedButtonRepresentation2D> CLOSEbuttonRepresentation;
	CLOSEbuttonRepresentation->SetNumberOfStates(1);
	CLOSEbuttonRepresentation->SetButtonTexture(0, Texture2);
	vtkNew<vtkButtonWidget> CLOSEbuttonWidget;
	CLOSEbuttonWidget->SetInteractor(renderWindowInteractor);
	CLOSEbuttonWidget->SetRepresentation(CLOSEbuttonRepresentation);
	double bds2[6] = {
		BTN2x - width, BTN2x + width,  // Adjust width boundaries
		BTN2y - height, BTN2y + height,  // Adjust height boundaries
		0.0, 0.0
	};
	CLOSEbuttonRepresentation->SetPlaceFactor(1);
	CLOSEbuttonRepresentation->PlaceWidget(bds2);
	CLOSEbuttonWidget->On();
	vtkNew<vtkTextActor> CLOSEbtntextActor;
	CLOSEbtntextActor->SetInput("EXIT");
	CLOSEbtntextActor->GetTextProperty()->SetFontSize(20);
	CLOSEbtntextActor->GetTextProperty()->SetBold(1);
	CLOSEbtntextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
	CenterTextActor(CLOSEbtntextActor, MainRenderer, BTN2x, BTN2y);
	MainRenderer->AddActor2D(CLOSEbtntextActor);
	vtkNew<vtkCallbackCommand> CLOSEbuttoncallback;
	CLOSEbuttoncallback->SetCallback(CLOSEButtonCallbackFunction);
	CLOSEbuttonWidget->AddObserver(vtkCommand::StateChangedEvent, CLOSEbuttoncallback);


	//// Add the mouse move callback
	//vtkNew<vtkCallbackCommand> mouseMoveCallback;
	//mouseMoveCallback->SetCallback(MouseMoveCallbackFunction);
	//renderWindowInteractor->AddObserver(vtkCommand::MouseMoveEvent, mouseMoveCallback);
	//// Setup the spline widget
	//vtkSmartPointer<vtkSplineWidget2> splineWidget = vtkSmartPointer<vtkSplineWidget2>::New();
	//splineWidget->SetInteractor(renderWindowInteractor);
	//splineWidget->CreateDefaultRepresentation();
	//splineWidget->On()


	// Camera setups
	vtkNew<vtkCameraOrientationWidget> camOrientManipulator;
	camOrientManipulator->SetParentRenderer(MainRenderer);
	camOrientManipulator->On();

	MainRenderer->ResetCamera();
	MainRenderer->GetActiveCamera()->SetPosition(0, 0, 160);
	MainRenderer->GetActiveCamera()->SetFocalPoint(0, 0, 0);
	MainRenderer->GetActiveCamera()->SetViewUp(0, 1, 0);
	MainRenderer->ResetCameraClippingRange();

	InsetRenderer->ResetCamera();
	InsetRenderer->GetActiveCamera()->SetPosition(95, 95, 40);
	InsetRenderer->GetActiveCamera()->SetFocalPoint(0, 0, 20);
	InsetRenderer->GetActiveCamera()->SetViewUp(0, 0, 1);
	InsetRenderer->ResetCameraClippingRange();

	// Custom Interaction Style
	vtkNew<C_InteractorStyle> style;
	style->SetDefaultRenderer(MainRenderer);
	style->MeshActor = movableActor;
	MainRenderer->AddActor(style->TextActor);
	renderWindowInteractor->SetInteractorStyle(style);

	renderWindow->Render();
	//renderWindowInteractor->ProcessEvents();

		// Get the HWND from the VTK Render Window
	HWND hWnd = static_cast<HWND>(renderWindow->GetGenericWindowId());
	SetWindowLong(hWnd, GWL_EXSTYLE, GetWindowLong(hWnd, GWL_EXSTYLE) | WS_EX_LAYERED); // | WS_EX_CLIENTEDGE
	SetLayeredWindowAttributes(hWnd, 0, WIN_TRANS, LWA_ALPHA);

	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();

	if (DEBUG) std::cout << Yellow << "        Viewer Exited." << ColorEnd << std::endl;

	Xoffset = style->X_offset;
	Yoffset = style->Y_offset;

	if (DEBUG) std::cout << Yellow << "        Offsets: " << ColorEnd << "X" << Xoffset << "  Y" << Yoffset
		<< "  Z" << CutHeight << "  Rot Z" << RotZ << std::endl;
	std::cout << std::endl;
}

void clean_difference(Cgal_Mesh& mesh, bool DEBUG) {
	auto vpm = get(CGAL::vertex_point, mesh);
	std::vector<std::size_t> component_ids(num_faces(mesh));
	auto component_map = CGAL::make_property_map(component_ids);
	std::size_t num_components = PMP::connected_components(mesh, component_map, PMP::parameters::vertex_point_map(vpm));

	std::vector<std::size_t> component_sizes(num_components, 0);
	for (Cgal_face_descriptor fd : faces(mesh)) {
		++component_sizes[component_map[fd]];
	}

	std::size_t largest_component_id = std::distance(component_sizes.begin(),
		std::max_element(component_sizes.begin(), component_sizes.end()));

	std::vector<Cgal_Face> faces_to_remove;
	for (Cgal_face_descriptor fd : faces(mesh)) {
		if (component_map[fd] != largest_component_id) {
			faces_to_remove.push_back(fd);
		}
	}

	for (Cgal_Face f : faces_to_remove) {
		mesh.remove_face(f);
	}

	mesh.collect_garbage();
	PMP::remove_isolated_vertices(mesh);
}

int close_mesh_hole(Cgal_Mesh& mesh, bool DEBUG) {
	int holes_closed = 0;
	std::vector<Cgal_Halfedge> border_halfedges;
	for (Cgal_Halfedge h : mesh.halfedges()) {
		if (mesh.is_border(h)) {
			border_halfedges.push_back(h);
		}
	}

	for (Cgal_Halfedge h : border_halfedges) {
		std::vector<Cgal_Face> patch_facets;
		std::vector<Cgal_Vertex> patch_vertices;

		bool success = std::get<0>(PMP::triangulate_refine_and_fair_hole(mesh,
			h,
			CGAL::parameters::face_output_iterator(std::back_inserter(patch_facets))
			.vertex_output_iterator(std::back_inserter(patch_vertices))));

		if (success) {
			holes_closed++;
		}
		else {
			std::cerr << Red << "        Failed to close mesh holes." << ColorEnd << std::endl;
		}
	}
	return holes_closed;
}

bool repair_and_validate_mesh(Cgal_Mesh& mesh, bool DEBUG) {
	bool isValid;
	Cgal_Mesh temp = mesh;
	mesh.clear();

	std::size_t Mesh_stitches = PMP::stitch_borders(temp);
	if (DEBUG) std::cout << Yellow << "        Mesh stitches: "
		<< ColorEnd << Mesh_stitches << std::endl;

	std::vector<Cgal_Face> degenerate_faces;
	PMP::degenerate_faces(temp, std::back_inserter(degenerate_faces));
	PMP::remove_degenerate_faces(degenerate_faces, temp);
	if (DEBUG) std::cout << Yellow << "        Removed degenerate faces: " << ColorEnd << degenerate_faces.size() << std::endl;

	std::size_t Mesh_New_vertices = PMP::duplicate_non_manifold_vertices(temp);
	if (DEBUG) std::cout << Yellow << "        Mesh new vertices: "
		<< ColorEnd << Mesh_New_vertices << std::endl;

	std::size_t Mesh_removed_vertices = PMP::remove_isolated_vertices(temp);
	if (DEBUG) std::cout << Yellow << "        Mesh removed vertices: "
		<< ColorEnd << Mesh_removed_vertices << std::endl;

	std::size_t removed_components = PMP::remove_connected_components_of_negligible_size(mesh);
	if (DEBUG) std::cout << Yellow << "        Removed small connected components: " << ColorEnd << removed_components << std::endl;

	PMP::orient(temp);

	if (!PMP::is_outward_oriented(temp)) {
		if (DEBUG) std::cout << Yellow << "        The mesh is not outward oriented." << ColorEnd << std::endl;
		PMP::reverse_face_orientations(temp); // Reverse the face orientations
		if (DEBUG) std::cout << Green << "        Orientation corrected." << ColorEnd << std::endl;
	}

	temp.collect_garbage();

	if (!CGAL::is_closed(temp))
		std::cerr << Red << "        Error: The mesh is not closed." << ColorEnd << std::endl;

	std::stringstream buffer;
	std::streambuf* prevcerr = std::cerr.rdbuf(buffer.rdbuf());
	isValid = CGAL::is_valid_polygon_mesh(temp, DEBUG);
	std::cerr.rdbuf(prevcerr);
	std::string line;
	while (std::getline(buffer, line)) {
		if (DEBUG) std::cout << Yellow << "        Validation output: " << ColorEnd << line << ColorEnd << std::endl;
	}

	mesh = temp;
	return isValid;
}

bool prepare_mesh(Cgal_Mesh& mesh) {
	if (!CGAL::is_triangle_mesh(mesh)) {
		std::cerr << "        Input mesh is not a triangle mesh." << std::endl;
	}
	assert(CGAL::is_valid_polygon_mesh(mesh));
	PMP::triangulate_faces(mesh);

	if (CGAL::is_closed(mesh)) printf("        Mesh: is closed\n");
	if (PMP::does_self_intersect(mesh)) printf("        Mesh: is self_intersect\n");

	PMP::remove_degenerate_faces(mesh);
	PMP::remove_isolated_vertices(mesh);
	PMP::duplicate_non_manifold_vertices(mesh);
	if (PMP::does_self_intersect(mesh)) {
		if (!PMP::experimental::remove_self_intersections(mesh,
			CGAL::parameters::preserve_genus(false))) {
			printf("        Mesh: Can not remove_self_intersections\n");
		}
	}

	//// close mesh
	//// fill hole
	//std::vector<Halfedge> border_cycles;
	//PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
	//for (Halfedge h : border_cycles) {
	//	std::vector<Face> patch_facets;
	//	PMP::triangulate_hole(
	//		mesh,
	//		h,
	//		std::back_inserter(patch_facets));
	//}

	if (CGAL::is_closed(mesh)) printf("        Mesh: is closed\n");
	else {
		std::cerr << "        Input mesh is not closed. Attempting to close the mesh." << std::endl;
		PMP::stitch_borders(mesh);
	}
	if (PMP::does_self_intersect(mesh)) printf("        Mesh: is self_intersect\n");

	//PMP::fair(mesh);
	std::cout << "        Smoothing the mesh..." << std::endl;
	std::set<Cgal_Mesh::Vertex_index> constrained_vertices;
	for (Cgal_Mesh::Vertex_index v : vertices(mesh)) {
		if (is_border(v, mesh))
			constrained_vertices.insert(v);
	}

	CGAL::Boolean_property_map<std::set<Cgal_Mesh::Vertex_index>> vcmap(constrained_vertices);
	PMP::smooth_shape(mesh, 0.0001, CGAL::parameters::number_of_iterations(10).vertex_is_constrained_map(vcmap));

	return true;
}

double get_height(Cgal_Mesh mesh) {
	std::vector<Cgal_Point> points;
	for (auto v : mesh.vertices()) {
		points.push_back(mesh.point(v));
	}
	Cgal_cuboid bbox = CGAL::bounding_box(points.begin(), points.end());
	double modelHeight = bbox.zmax() - bbox.zmin();
	if (DEBUG) std::cout << Yellow << "        Mesh Height:  H " << ColorEnd << modelHeight << std::endl;
	return modelHeight;
}

void get_mesh_dimensions(Cgal_Mesh mesh, double& modelWidth, double& modelLength, double& modelHeight, bool DEBUG) {
	std::vector<Cgal_Point> points;
	for (auto v : mesh.vertices()) {
		points.push_back(mesh.point(v));
	}
	Cgal_cuboid bbox = CGAL::bounding_box(points.begin(), points.end());
	modelWidth = bbox.xmax() - bbox.xmin();
	modelLength = bbox.ymax() - bbox.ymin();
	modelHeight = bbox.zmax() - bbox.zmin();
	if (DEBUG) std::cout << Yellow << "        Mesh Dimensions:  W " << ColorEnd
		<< modelWidth << "   L "
		<< modelLength << "   H "
		<< modelHeight << std::endl;
}

void get_mesh_center(Cgal_Mesh mesh, Cgal_Point& center, bool DEBUG) {
	CGAL::Bbox_3 bbox;
	for (auto v : mesh.vertices()) {
		bbox += mesh.point(v).bbox();
	}
	center = Cgal_Point((bbox.xmin() + bbox.xmax()) / 2.0, (bbox.ymin() + bbox.ymax()) / 2.0, (bbox.zmin() + bbox.zmax()) / 2.0);
	if (DEBUG) std::cout << Yellow << "        Mesh Center:  " << ColorEnd
		<< center.x() << "  "
		<< center.y() << "  "
		<< center.z() << std::endl;
}

void get_mesh_centroid(Cgal_Mesh mesh, Cgal_Point& centroid, bool DEBUG) {
	std::vector<Cgal_Point> vertices;
	for (auto v : mesh.vertices()) {
		vertices.push_back(mesh.point(v));
	}
	centroid = CGAL::centroid(vertices.begin(), vertices.end());
	if (DEBUG) std::cout << Yellow << "        Mesh Centroid:  " << ColorEnd
		<< centroid.x() << "  "
		<< centroid.y() << "  "
		<< centroid.z() << std::endl;
}

void settle_mesh(Cgal_Mesh& mesh, bool DEBUG) {
	double min_z = std::numeric_limits<double>::infinity();
	for (auto v : mesh.vertices()) {
		double z = mesh.point(v).z();
		if (z < min_z) min_z = z;
	}
	Cgal_Vector translation_vector(0, 0, -min_z);
	for (auto v : mesh.vertices()) {
		Cgal_Point p = mesh.point(v) + translation_vector;
		mesh.point(v) = p;
	}
	if (DEBUG) std::cout << Yellow << "        Mesh Settled at Z:  " << ColorEnd << -min_z << std::endl;
}

void cut_mesh(Cgal_Mesh& mesh, double height, bool DEBUG) {
	if (height <= -0.0001) {
		Cgal_Plane plane(0, 0, -1, -height); // Adjust 'distance' as needed
		if (PMP::clip(mesh, plane, PMP::parameters::clip_volume(true))) {
			double min_z = std::numeric_limits<double>::infinity();
			for (auto v : mesh.vertices()) {
				double z = mesh.point(v).z();
				if (z < min_z) min_z = z;
			}
			Cgal_Vector translation_vector(0, 0, min_z);
			for (auto v : mesh.vertices()) {
				Cgal_Point p = mesh.point(v) - translation_vector;
				mesh.point(v) = p;
			}
			if (DEBUG) std::cout << Yellow << "        Mesh Cut at Z:  " << ColorEnd << height << "  , Settled at Z: " << min_z << std::endl;
		}
		else {
			std::cerr << Red << "        Cutting mesh failed." << ColorEnd << std::endl;
			if (DEBUG) std::cout << Yellow << "        Trying another Cutting way..." << ColorEnd << std::endl;
		
			double size = 100.0, bottom_z = -20;
			Cgal_Mesh clipper, Result_Mesh;
			Cgal_Vertex v0 = clipper.add_vertex(Cgal_Point(-size, -size, -height));
			Cgal_Vertex v1 = clipper.add_vertex(Cgal_Point(size, -size, -height));
			Cgal_Vertex v2 = clipper.add_vertex(Cgal_Point(size, size, -height));
			Cgal_Vertex v3 = clipper.add_vertex(Cgal_Point(-size, size, -height));
			Cgal_Vertex v4 = clipper.add_vertex(Cgal_Point(-size, -size, bottom_z));
			Cgal_Vertex v5 = clipper.add_vertex(Cgal_Point(size, -size, bottom_z));
			Cgal_Vertex v6 = clipper.add_vertex(Cgal_Point(size, size, bottom_z));
			Cgal_Vertex v7 = clipper.add_vertex(Cgal_Point(-size, size, bottom_z));
			clipper.add_face(v0, v1, v2); clipper.add_face(v2, v3, v0); // Top face
			clipper.add_face(v4, v6, v5); clipper.add_face(v6, v4, v7); // Bottom face
			clipper.add_face(v0, v4, v1); clipper.add_face(v1, v4, v5); // Four side faces
			clipper.add_face(v1, v5, v2); clipper.add_face(v2, v5, v6);
			clipper.add_face(v2, v6, v3); clipper.add_face(v3, v6, v7);
			clipper.add_face(v3, v7, v0); clipper.add_face(v0, v7, v4);
		
			if (PMP::corefine_and_compute_difference(mesh, clipper, Result_Mesh)) {
				double min_z = std::numeric_limits<double>::infinity();
				for (auto v : mesh.vertices()) {
					double z = mesh.point(v).z();
					if (z < min_z) min_z = z;
				}
				Cgal_Vector translation_vector(0, 0, min_z);
				for (auto v : mesh.vertices()) {
					Cgal_Point p = mesh.point(v) - translation_vector;
					mesh.point(v) = p;
				}
				if (DEBUG) std::cout << Yellow << "        Mesh Cut at Z:  " << ColorEnd << height << "  , Settled at Z: " << min_z << std::endl;
				mesh.clear();
				mesh = Result_Mesh;
			}
			else {
				std::cerr << Red << "        Cutting mesh failed." << ColorEnd << std::endl;
			}
		}
	}
}

void extrude_mesh(Cgal_Mesh& mesh, double target_z, bool DEBUG) {
	if (target_z >= 0.0001) {
		double min_z = 0.001;
		for (Cgal_Vertex v : mesh.vertices()) {
			Cgal_Point& p = mesh.point(v);
			if (p.z() >= min_z) mesh.point(v) = Cgal_Point(p.x(), p.y(), p.z() + target_z);
		}
		if (DEBUG) std::cout << Yellow << "        Mesh Extruded:  " << ColorEnd << target_z << std::endl;
		settle_mesh(mesh, DEBUG);
	}
}

void scale_mesh(Cgal_Mesh& mesh, double XYscale, double XYtopscale, double Zscale, double zThreshold, bool DEBUG) {
	for (auto v : mesh.vertices()) {
		Cgal_Point& point = mesh.point(v);
		double new_x, new_y, new_z;
		if (point.z() > zThreshold) {
			new_x = point.x() * XYtopscale;
			new_y = point.y() * XYtopscale;
		}
		else {
			new_x = point.x() * XYscale;
			new_y = point.y() * XYscale;
		}
		new_z = point.z() * Zscale;
		mesh.point(v) = Cgal_Point(new_x, new_y, new_z);
	}
}

void translate_mesh(Cgal_Mesh& mesh, double Vx, double Vy, double Vz, bool DEBUG) {
	Cgal_Transformation translation(CGAL::TRANSLATION, Cgal_Vector(Vx, Vy, Vz));
	for (auto v : mesh.vertices()) {
		mesh.point(v) = translation(mesh.point(v));
	}
	if (DEBUG) std::cout << Yellow << "        translation Applied:  " << ColorEnd << Cgal_Vector(Vx, Vy, Vz) << std::endl;
}

void rotate_mesh(Cgal_Mesh& mesh, double x_deg, double y_deg, double z_deg, bool DEBUG) {
	if (z_deg < -0.001 || z_deg > 0.001) {
		double cos_x = std::cos(x_deg * M_PI / 180.0), sin_x = std::sin(x_deg * M_PI / 180.0);
		double cos_y = std::cos(y_deg * M_PI / 180.0), sin_y = std::sin(y_deg * M_PI / 180.0);
		double cos_z = std::cos(z_deg * M_PI / 180.0), sin_z = std::sin(z_deg * M_PI / 180.0);

		Cgal_Transformation rot_mtx_x(1, 0, 0, 0, 0, cos_x, -sin_x, 0, 0, sin_x, cos_x, 0, 1);
		Cgal_Transformation rot_mtx_y(cos_y, 0, sin_y, 0, 0, 1, 0, 0, -sin_y, 0, cos_y, 0, 1);
		Cgal_Transformation rot_mtx_z(cos_z, -sin_z, 0, 0, sin_z, cos_z, 0, 0, 0, 0, 1, 0, 1);
		Cgal_Transformation combined = rot_mtx_z * rot_mtx_y * rot_mtx_x;
		PMP::transform(combined, mesh);
		if (DEBUG) std::cout << Yellow << "        Rotation Applied:  "
			<< ColorEnd << "X " << x_deg << ", Y " << y_deg << ", Z " << z_deg << std::endl;
	}
}

void translate_pts(std::vector<Cgal_Point>& points, double Vx, double Vy, double Vz, bool DEBUG) {
	Cgal_Transformation translation(CGAL::TRANSLATION, Cgal_Vector(Vx, Vy, Vz));
	for (auto& point : points) {
		point = translation(point);
	}
	if (DEBUG) std::cout << Yellow << "        PTS translation Applied:  " << ColorEnd << Cgal_Vector(Vx, Vy, Vz) << std::endl;
}

void rotate_pts(std::vector<Cgal_Point>& points, double x_deg, double y_deg, double z_deg, bool DEBUG) {
	if (z_deg < -0.001 || z_deg > 0.001) {
		double cos_x = std::cos(x_deg * M_PI / 180.0), sin_x = std::sin(x_deg * M_PI / 180.0);
		double cos_y = std::cos(y_deg * M_PI / 180.0), sin_y = std::sin(y_deg * M_PI / 180.0);
		double cos_z = std::cos(z_deg * M_PI / 180.0), sin_z = std::sin(z_deg * M_PI / 180.0);

		Cgal_Transformation rot_mtx_x(1, 0, 0, 0, 0, cos_x, -sin_x, 0, 0, sin_x, cos_x, 0, 1);
		Cgal_Transformation rot_mtx_y(cos_y, 0, sin_y, 0, 0, 1, 0, 0, -sin_y, 0, cos_y, 0, 1);
		Cgal_Transformation rot_mtx_z(cos_z, -sin_z, 0, 0, sin_z, cos_z, 0, 0, 0, 0, 1, 0, 1);
		Cgal_Transformation combined = rot_mtx_z * rot_mtx_y * rot_mtx_x;
		for (auto& point : points) {
			point = combined(point);
		}
		if (DEBUG) std::cout << Yellow << "        PTS Rotation Applied:  "
			<< ColorEnd << "X " << x_deg << ", Y " << y_deg << ", Z " << z_deg << std::endl;
	}
}

bool read_STL_data(std::string identifier, Cgal_Mesh& mesh, bool DEBUG) {
	auto start_ = std::chrono::high_resolution_clock::now();
	mesh.clear();
	for (const auto& data : FONT_STL) {
		if (data.key == identifier) {
			std::istringstream iss(std::string(reinterpret_cast<const char*>(data.data), data.size), std::ios::binary);
			if (CGAL::IO::read_STL(iss, mesh)) {
				auto finish_ = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed_ = finish_ - start_;
				if (DEBUG) std::cout << "     >> " << Yellow << "Read Mesh: " << ColorEnd << identifier
					<< Cyan << "   ET: " << elapsed_.count() << " Sec" << ColorEnd << std::endl;
				return true;
			}
			break;
		}
	}
	std::cerr << Red << "        Error: No STL data available for:  " << ColorEnd << identifier << std::endl;
	return false;
}

bool read_STL_Cgal(std::string filename, Cgal_Mesh& mesh, bool DEBUG) {
	auto start_ = std::chrono::high_resolution_clock::now();
	mesh.clear();
	fs::path filepath(filename);
	////if (!CGAL::IO::read_STL(filename, mesh)) {
	//if (!PMP::IO::read_polygon_mesh(filename, mesh)) {
	//	std::cerr << Red << "        Error: Cannot read the STL file:  " << ColorEnd << filepath.filename().string() << std::endl;
	//	return false;
	//}

	std::vector<Cgal_Point> points;
	std::vector<std::vector<std::size_t>> polygons;
	if (!CGAL::IO::read_polygon_soup(filename, points, polygons)) return false;
	PMP::repair_polygon_soup(points, polygons);
	PMP::orient_polygon_soup(points, polygons);
	PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);

	auto finish_ = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_ = finish_ - start_;
	if (DEBUG) std::cout << "     >> " << Yellow << "Read STL File: " << ColorEnd << filepath.filename().string()
		<< Cyan << "   ET: " << elapsed_.count() << " Sec" << ColorEnd << std::endl;
	return true;
}

bool write_STL_Cgal(std::string filename, Cgal_Mesh mesh, bool DEBUG) {
	auto start_ = std::chrono::high_resolution_clock::now();
	fs::path filepath(filename);
	//if (!CGAL::IO::write_polygon_mesh(filename, mesh, CGAL::parameters::stream_precision(10))) {
	if (!CGAL::IO::write_polygon_mesh(filename, mesh)) {
		std::cerr << Red << "        Error: Cannot write the STL file:  " << ColorEnd << filepath.filename().string() << std::endl;
		return false;
	}
	auto finish_ = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_ = finish_ - start_;
	if (DEBUG) std::cout << "     << " << Yellow << "Written STL File: " << ColorEnd << filepath.filename().string()
		<< Cyan << "   ET: " << elapsed_.count() << " Sec" << ColorEnd << std::endl;
	return true;
}

void remesh_mesh(Cgal_Mesh& mesh, double target_edge_length, unsigned int number_of_iterations, bool DEBUG) {
	auto start_ = std::chrono::high_resolution_clock::now();
	if (DEBUG) std::cout << Yellow << "        Creating a remesh" << ColorEnd << std::endl;
	PMP::triangulate_faces(mesh);
	if (DEBUG) std::cout << Yellow << "        Triangulate faces" << ColorEnd << std::endl;
	PMP::remove_degenerate_faces(mesh);
	if (DEBUG) std::cout << Yellow << "        Remove degenerate faces" << ColorEnd << std::endl;
	PMP::stitch_borders(mesh);
	if (DEBUG) std::cout << Yellow << "        Stitch borders" << ColorEnd << std::endl;

	std::vector<std::pair<Cgal_face_descriptor, Cgal_face_descriptor>> self_intersections;
	PMP::self_intersections(mesh, std::back_inserter(self_intersections));
	if (!self_intersections.empty()) {
		PMP::experimental::autorefine_and_remove_self_intersections(mesh);
		if (DEBUG) std::cout << Yellow << "        Autorefine and remove self intersections" << ColorEnd << std::endl;
	}

#pragma warning(push)
#pragma warning(disable: 4996) 
	bool print = false;
	for (auto h : mesh.halfedges()) {
		if (mesh.is_border(h)) {
			std::vector<Cgal_face_descriptor> patch;
			PMP::triangulate_hole(mesh, h, std::back_inserter(patch));
			if (DEBUG && !print) std::cout << Yellow << "        Triangulate hole" << ColorEnd << std::endl;
			print = true;
		}
	}
#pragma warning(pop)

	//PMP::isotropic_remeshing(mesh.faces(),
	//	target_edge_length,
	//	mesh,
	//	PMP::parameters::number_of_iterations(number_of_iterations));

	auto finish_ = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_ = finish_ - start_;
	if (DEBUG) std::cout << "        " << Yellow << "Remesh done " << ColorEnd
		<< Cyan << "   ET: " << elapsed_.count() << " Sec" << ColorEnd << std::endl;
}

void create_fixture(std::string ID_Str, Cgal_Mesh Fix_Mesh, Cgal_Mesh& Res_Mesh, bool DEBUG) {
	bool lastWasDigit = false;
	double offsetX = -6.5, offsetY = -7.5, offsetZ = 4.0;
	double XYscale = 0.18, XYtopscale = 0.18, Zscale = 0.30;
	double zThreshold = 0.1;
	double Xspacing = 0.8, Yspacing = 2.9;
	double zDepth = -0.7;
	Cgal_Mesh Tag_Mesh;

	if (DEBUG) std::cout << Yellow << "        Creating tag on Fixture:  " << ColorEnd << ID_Str << std::endl;

	std::transform(ID_Str.begin(), ID_Str.end(), ID_Str.begin(), [](unsigned char c) { return std::toupper(c); });
	for (char c : ID_Str) {
		Cgal_Mesh Letter_Mesh;
		double FontWidth = 0.0, FontLength = 0.0, FontHeight = 0.0;

		if (!read_STL_data(std::string(1, c), Letter_Mesh, false)) continue;

		get_mesh_dimensions(Letter_Mesh, FontWidth, FontLength, FontHeight, false);

		if (std::isdigit(c)) lastWasDigit = true;
		else if (lastWasDigit) {
			offsetY -= (FontLength * XYscale) + Yspacing;
			offsetX = -6.35; // 0.15
			lastWasDigit = false;
		}

		scale_mesh(Letter_Mesh, XYscale, XYtopscale, Zscale, zThreshold, false);
		translate_mesh(Letter_Mesh, offsetX, offsetY, offsetZ + zDepth, false);
		offsetX += (FontWidth * XYscale) + Xspacing;
		CGAL::copy_face_graph(Letter_Mesh, Tag_Mesh);
	}

	Res_Mesh.clear();
	if (!PMP::corefine_and_compute_difference(Fix_Mesh, Tag_Mesh, Res_Mesh))
		std::cerr << Red << "        Tag Subtraction failed for: " << ColorEnd << ID_Str << std::endl;
}


bool modify_model_mesh(Cgal_Mesh F_Mesh, std::string I_Path, std::string O_Path, const std::string ID_,
	double M_Xoffset, double M_Yoffset, double C_height, double M_Zrot, bool DEBUG) {
	auto start_ = std::chrono::high_resolution_clock::now();

	Cgal_Mesh M_Mesh, F_ID_Mesh, Result_Mesh;
	if (DEBUG) std::cout << "      > " << Yellow << "Applied Offsets: " << ColorEnd
		<< "X" << M_Xoffset << "  Y" << M_Yoffset << "  CH" << C_height << "  R" << M_Zrot << std::endl;

	read_STL_Cgal(I_Path, M_Mesh, DEBUG);

	remesh_mesh(M_Mesh, 0.3, 3, DEBUG);
	
	create_fixture(ID_, F_Mesh, F_ID_Mesh, DEBUG);

	rotate_mesh(M_Mesh, 0, 0, M_Zrot, DEBUG);

	Cgal_Point M_center;
	get_mesh_center(M_Mesh, M_center, DEBUG);
	translate_mesh(M_Mesh, M_Xoffset - M_center.x(), M_Yoffset - M_center.y(), 0, DEBUG);

	extrude_mesh(M_Mesh, C_height, DEBUG);
	cut_mesh(M_Mesh, C_height, DEBUG);

	//CGAL::copy_face_graph(F_ID_Mesh, Result_Mesh);
	//CGAL::copy_face_graph(M_Mesh, Result_Mesh);
	//if (DEBUG) std::cout << Yellow << "        Fixture added to: " << ColorEnd << ID_ << std::endl;

	PMP::corefine_and_compute_union(M_Mesh, F_ID_Mesh, Result_Mesh);
	if (DEBUG) std::cout << Yellow << "        Fixture added to: " << ColorEnd << ID_ << std::endl;


	//prepare_mesh(Result_Mesh);
	//repair_and_validate_mesh(Result_Mesh, DEBUG);


	if (!write_STL_Cgal(O_Path, Result_Mesh, DEBUG)) return false;

	auto finish_ = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_ = finish_ - start_;
	std::cout << "   >><< " << Green << "Operation completed successfully." << ColorEnd
		<< Cyan << "   ET: " << elapsed_.count() << " Sec" << ColorEnd << std::endl;
	std::cout << std::endl;
	return true;
}

void process_files(std::vector<std::string> filenames, std::string filter, std::string inputPath,
	std::string outputPath, bool DEBUG) {

	std::cout << std::endl;
	std::vector<std::string> filtered_files;
	std::copy_if(filenames.begin(), filenames.end(), std::back_inserter(filtered_files),
		[&filter](const std::string& name) {
			return name.find(filter) != std::string::npos;
		});
	std::sort(filtered_files.begin(), filtered_files.end());

	Cgal_Mesh Fixture_Mesh;
	read_STL_data("fixture", Fixture_Mesh, DEBUG);

	Cgal_Point V_center;
	double Model_Xoffset = 0.0, Model_Yoffset = 0.0, cut_height = 0.0, Model_Zrot = 0.0;
	bool started = false;
	for (const auto& fileID : filtered_files) {
		std::string Model_In_Path = inputPath + fileID + ".stl";
		std::string Model_Out_Path = outputPath + fileID + ".stl";

		model_count++;
		std::cout << (model_count < 10 ? "   0" : "   ") + std::to_string(model_count)
			<< " > " << Yellow << "Modifying " << ColorEnd << fileID << std::endl;

		if (fileID == filtered_files[0]) {
			Cgal_Mesh Visual_mesh;
			read_STL_Cgal(Model_In_Path, Visual_mesh, DEBUG);
			get_mesh_center(Visual_mesh, V_center, DEBUG);
			translate_mesh(Visual_mesh, -V_center.x(), -V_center.y(), 0, DEBUG);
			visualize_mesh(Visual_mesh, get_height(Visual_mesh), fileID, Fixture_Mesh,
				Model_Xoffset, Model_Yoffset, cut_height, Model_Zrot, true);
		}
		modify_model_mesh(Fixture_Mesh, Model_In_Path, Model_Out_Path, fileID,
			Model_Xoffset, Model_Yoffset, cut_height, Model_Zrot, DEBUG);
	}
}

void start_viewer_process(std::map<std::string, std::vector<std::string>> groupedFiles, std::vector<std::string> sortedKeys, std::string& inputPath, std::string& outputPath) {
	for (const auto& key : sortedKeys) {
		std::cout << Green << "        Processing CaseID: " << ColorEnd << key << std::endl;
		int Count = 0;
		for (const auto& filename : groupedFiles[key]) {
			Count++;
			std::cout << (Count < 10 ? "   0" : "   ") + std::to_string(Count) << " > " << Cyan << filename << ColorEnd << std::endl;
		}

		// LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER
		process_files(groupedFiles[key], "L", inputPath, outputPath, DEBUG);
		// UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER 
		process_files(groupedFiles[key], "U", inputPath, outputPath, DEBUG);
	}
}


int main(int argc, char* argv[]) {
	setConsoleSize(73, 35);

	std::cout << Cyan << "\n===========================" << ColorEnd
		<< Yellow << "'Created by Banna'" << ColorEnd
		<< Cyan << "===========================" << ColorEnd << std::endl;

	std::cout << Cyan << "======================" << ColorEnd
		<< Yellow << "'AB FIXTURE CREATOR TOOL V3'" << ColorEnd
		<< Cyan << "======================" << ColorEnd << std::endl;
	std::cout << Cyan << "========================================================================\n" << ColorEnd << std::endl;

	auto start = std::chrono::high_resolution_clock::now();
	std::string Main_path = fs::current_path().string();
	std::string outputPath = Main_path + "/output/";
	std::string inputPath = Main_path + "/input/";
	bool has_stl_files = false;
	std::vector<std::string> filenames;

	if (!fs::exists(inputPath)) fs::create_directory(inputPath);
	for (const auto& entry : fs::directory_iterator(inputPath)) {
		if (entry.path().extension() == ".stl") {
			std::string filename_WE = entry.path().stem().string();
			filenames.push_back(filename_WE);
			has_stl_files = true;
		}
	}

	if (!has_stl_files) {
		std::cout << Red << "        No STL files in input folder\n" << ColorEnd << std::endl;
		std::cout << "        Press enter to continue...";
		std::cin.get();
		return EXIT_SUCCESS;
	}

	if (!fs::exists(outputPath))
		fs::create_directory(outputPath);
	else for (const auto& entry : fs::directory_iterator(outputPath))
		fs::remove_all(entry.path());


	std::regex filenameRegex("([0-9]+)([a-zA-Z]+)([0-9]+)");
	std::map<std::string, std::vector<std::string>> groupedFiles;

	for (const auto& filename : filenames) {
		std::smatch matches;
		if (std::regex_match(filename, matches, filenameRegex) && matches.size() > 1) {
			std::string groupKey = matches[1].str();  // First numeric sequence
			groupedFiles[groupKey].push_back(filename);
		}
	}


	std::vector<std::string> sortedKeys;
	for (const auto& pair : groupedFiles)
		sortedKeys.push_back(pair.first);
	//std::sort(sortedKeys.begin(), sortedKeys.end());
	std::sort(sortedKeys.begin(), sortedKeys.end(), [](const std::string& a, const std::string& b) {
		return std::stoi(a) < std::stoi(b);
		});


	start_viewer_process(groupedFiles, sortedKeys, inputPath, outputPath);

	//for (const auto& key : sortedKeys) {
	//	std::cout << Green << "        Processing CaseID: " << ColorEnd << key << std::endl;
	//	int Count = 0;
	//	for (const auto& filename : groupedFiles[key]) {
	//		Count++;
	//		std::cout << (Count < 10 ? "   0" : "   ") + std::to_string(Count) << " > " << Cyan << filename << ColorEnd << std::endl;
	//	}
	//
	//	// LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER LOWER
	//	process_files(groupedFiles[key], "L", inputPath, outputPath, DEBUG);
	//	// UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER UPPER 
	//	process_files(groupedFiles[key], "U", inputPath, outputPath, DEBUG);
	//}

	std::cout << std::endl;
	displayUserName();

	auto finish = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(finish - start);
	int minutes = elapsed.count() / 60;
	int seconds = elapsed.count() % 60;
	std::cout << Yellow << "        Elapsed time: " << minutes << " Minutes " << seconds << " Seconds" << ColorEnd << std::endl;

	std::cout << std::endl;
	std::cout << "        Press enter to continue...";
	std::cin.get();  // Waits for the user to press Enter
	return EXIT_SUCCESS;
}