#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <QDebug>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <vector>
#include <cmath>
#include <limits>


namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    EdgeAttributes( OpenMesh::Attributes::Color | OpenMesh::Attributes::Status );
    // vertex thickness
    VertexTraits{float thickness; float value;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

struct VertexMatching{
    MyMesh::VertexIter v1;
    MyMesh::VertexIter v2;
};

typedef VertexMatching VMatch;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void Match();
    void displayMesh(MyMesh *_mesh, int widget_nbr, bool isTemperatureMap = false, float mapRange = -1);
    void resetAllColorsAndThickness(MyMesh* _mesh);

private slots:
    void on_pushButton_chargement_clicked();
    void on_pushButton_chargement2_clicked();

    void on_pushButton_morphing_clicked();

    void on_animationBar_sliderReleased();

    void on_button_parametrization_clicked();

    void on_button_parametrization2_clicked();

private:

    MyMesh mesh;
    MyMesh sphere1;
    OpenMesh::Vec3f center_m1;
    float max_m1;

    MyMesh mesh2;
    MyMesh sphere2;
    OpenMesh::Vec3f center_m2;
    float max_m2;

    MyMesh result;

    std::vector<VMatch> transformation;
    float square(float a);
    float dist(OpenMesh::Vec3f p1, OpenMesh::Vec3f p2);
    void Normalization_mesh1();
    void Normalization_mesh2();
    void Normalizations();
    void on_pushButton_generer_clicked();



    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
