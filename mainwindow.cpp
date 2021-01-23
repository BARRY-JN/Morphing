#include "mainwindow.h"
#include "ui_mainwindow.h"

float MainWindow::square(float a){return a*a;}

float MainWindow::dist(OpenMesh::Vec3f p1, OpenMesh::Vec3f p2){
     //sqrt[(Xa-Xb)²+(Ya-Yb)²+(Za-Zb)²]
    return sqrtf(square(p1[0]-p2[0])+square(p1[1]-p2[1])+square(p1[2]-p2[2]));
}

void MainWindow::Normalization_mesh1(){
    OpenMesh::Vec3f p1;
    float mx=0,my=0,mz=0, total=0;

    //calcul du centre du modele
    for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it){
        p1 = mesh.point(*v_it);
        mx+=p1[0];
        my+=p1[1];
        mz+=p1[2];
    }
    mx=mx/mesh.n_vertices();
    my=my/mesh.n_vertices();
    mz=mz/mesh.n_vertices();

    //trouver le point le plus éloigné du centre
    center_m1[0]=mx;
    center_m1[1]=my;
    center_m1[2]=mz;


    max_m1=0;
    for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it){
        p1 = mesh.point(*v_it);

        float d = dist(center_m1,p1);
        if(d>max_m1)
            max_m1=d;

    }

    float radius=1;
    //parametrization
    for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it){
        p1 = mesh.point(*v_it);

        //vecteur deplacement
        float dx = p1[0]-center_m1[0];
        float dy = p1[1]-center_m1[1];
        float dz = p1[2]-center_m1[2];

        //vecteur unitaire
        float rx=(dx/max_m1);
        float ry=(dy/max_m1);
        float rz=(dz/max_m1);

        OpenMesh::Vec3f interpolated(rx,ry,rz);
        mesh.set_point( v_it,interpolated);
    }

    mesh.update_normals();
    resetAllColorsAndThickness(&mesh);
    displayMesh(&mesh, 1);
}

void MainWindow::Normalization_mesh2(){
    OpenMesh::Vec3f p1;
    float mx=0,my=0,mz=0, total=0;

    //calcul du centre du modele
    for (MyMesh::VertexIter v_it=mesh2.vertices_begin(); v_it!=mesh2.vertices_end(); ++v_it){
        p1 = mesh2.point(*v_it);
        mx+=p1[0];
        my+=p1[1];
        mz+=p1[2];
    }
    mx=mx/mesh2.n_vertices();
    my=my/mesh2.n_vertices();
    mz=mz/mesh2.n_vertices();

    //trouver le point le plus éloigné du centre
    center_m2[0]=mx;
    center_m2[1]=my;
    center_m2[2]=mz;


    max_m2=0;
    for (MyMesh::VertexIter v_it=mesh2.vertices_begin(); v_it!=mesh2.vertices_end(); ++v_it){
        p1 = mesh2.point(*v_it);

        float d = dist(center_m2,p1);
        if(d>max_m2)
            max_m2=d;

    }

    float radius=1;
    //parametrization
    for (MyMesh::VertexIter v_it=mesh2.vertices_begin(); v_it!=mesh2.vertices_end(); ++v_it){
        p1 = mesh2.point(*v_it);

        //vecteur deplacement
        float dx = p1[0]-center_m2[0];
        float dy = p1[1]-center_m2[1];
        float dz = p1[2]-center_m2[2];

        //vecteur unitaire
        float rx=(dx/max_m2);
        float ry=(dy/max_m2);
        float rz=(dz/max_m2);

        OpenMesh::Vec3f interpolated(rx,ry,rz);
        mesh2.set_point( v_it,interpolated);
    }

    mesh2.update_normals();
    resetAllColorsAndThickness(&mesh2);
    displayMesh(&mesh2, 2);
}
void MainWindow::Normalizations(){
    Normalization_mesh1();
    Normalization_mesh2();
}

void MainWindow::Match(){
     OpenMesh::Vec3f pointA, pointB;
     float min; //minimum=la plus grande valeur stockable dans un float
     transformation.clear();
     MyMesh::VertexIter saveIterator; //iterator qui contient la référence du point le plus proche actuel

     float i=0;
    for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it){
        pointA = mesh.point(*v_it);
        min=std::numeric_limits<float>::max();
        for (MyMesh::VertexIter v_it2=mesh2.vertices_begin(); v_it2!=mesh2.vertices_end(); ++v_it2){

            pointB = mesh2.point(*v_it2);
            float d=dist(pointA,pointB);

            if(d<min){
                min=d;
                saveIterator=v_it2;
            }
        }
        VMatch v; v.v1=v_it; v.v2=saveIterator;
        transformation.push_back(v);
        i++;
        ui->progressBar->setValue((int)((i/mesh.n_vertices())*100.0f));
    }
}

void MainWindow::on_pushButton_morphing_clicked()
{
    OpenMesh::Vec3f p1,p2;
    ui->pushButton_morphing->setEnabled(false);
    ui->animationBar->setEnabled(false);
    Match();
    ui->pushButton_morphing->setEnabled(true);
    ui->animationBar->setEnabled(true);





}


void MainWindow::on_animationBar_sliderReleased()
{
    OpenMesh::Vec3f p1,p2;
    float y1,y2,x1,x2,z1,z2,rx,ry,rz,pas;
    float value = ui->animationBar->value()/100.0;

    for(int i=0;i<transformation.size();i++){
        p1=mesh.point(*(transformation[i].v1)); x1=p1[0];y1=p1[1];z1=p1[2];
        p2=mesh2.point(*(transformation[i].v2)); x2=p2[0];y2=p2[1];z2=p2[2];

        //on effectue l'interpolation entre le double tableau qui contient l'association des points entre les 2 meshs
        float dx = x2-x1;
        float dy = y2-y1;
        float dz = z2-z1;

        rx=x1+dx*value;
        ry=y1+dy*value;
        rz=z1+dz*value;

        OpenMesh::Vec3f interpolated(rx,ry,rz);

        result.set_point( *(transformation[i].v1),interpolated);
    }
    result.update_normals();
    resetAllColorsAndThickness(&result);
    displayMesh(&result, 3);
}


void MainWindow::on_button_parametrization_clicked()
{
    OpenMesh::Vec3f p1;
    float radius=1;
    //parametrization
    for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it){
        p1 = mesh.point(*v_it);

        //vecteur deplacement
        float dx = p1[0]-center_m1[0];
        float dy = p1[1]-center_m1[1];
        float dz = p1[2]-center_m1[2];

        //on calcule combien de fois multiplier le vecteur unitaire pour avoir le rayon voulu
        float multi = radius/sqrtf(square(dx)+square(dy)+square(dz));

        float rx=dx*multi;
        float ry=dy*multi;
        float rz=dz*multi;

        OpenMesh::Vec3f interpolated(rx,ry,rz);
        mesh.set_point( v_it,interpolated);
    }

    mesh.update_normals();
    resetAllColorsAndThickness(&mesh);
    displayMesh(&mesh, 1);
}

void MainWindow::on_button_parametrization2_clicked()
{
    OpenMesh::Vec3f p1;

    float radius=1;
    //parametrization
    for (MyMesh::VertexIter v_it=mesh2.vertices_begin(); v_it!=mesh2.vertices_end(); ++v_it){
        p1 = mesh2.point(*v_it);

        //vecteur deplacement
        float dx = p1[0]-center_m2[0];
        float dy = p1[1]-center_m2[1];
        float dz = p1[2]-center_m2[2];

        //on calcule combien de fois multiplier le vecteur unitaire pour avoir le rayon voulu
        float multi = radius/sqrtf(square(dx)+square(dy)+square(dz));

        float rx=dx*multi;
        float ry=dy*multi;
        float rz=dz*multi;

        OpenMesh::Vec3f interpolated(rx,ry,rz);
        mesh2.set_point( v_it,interpolated);
    }

    mesh2.update_normals();
    resetAllColorsAndThickness(&mesh2);
    displayMesh(&mesh2, 2);
}
/* **** début de la partie boutons et IHM **** */


// exemple pour charger un fichier .obj
void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());
    OpenMesh::IO::read_mesh(result, fileName.toUtf8().constData());

    Normalization_mesh1();
    sphere1 = MyMesh(mesh);

    result.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&result);

    // on affiche le maillage
    displayMesh(&result, 3);
}

void MainWindow::on_pushButton_chargement2_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));
    OpenMesh::IO::read_mesh(mesh2, fileName.toUtf8().constData());
    Normalization_mesh2();
    sphere2 = MyMesh(mesh2);
}

// exemple pour construire un mesh face par face
void MainWindow::on_pushButton_generer_clicked()
{
    MyMesh mesh;

    // on construit une liste de sommets
    MyMesh::VertexHandle sommets[8];
    sommets[0] = mesh.add_vertex(MyMesh::Point(-1, -1,  1));
    sommets[1] = mesh.add_vertex(MyMesh::Point( 1, -1,  1));
    sommets[2] = mesh.add_vertex(MyMesh::Point( 1,  1,  1));
    sommets[3] = mesh.add_vertex(MyMesh::Point(-1,  1,  1));
    sommets[4] = mesh.add_vertex(MyMesh::Point(-1, -1, -1));
    sommets[5] = mesh.add_vertex(MyMesh::Point( 1, -1, -1));
    sommets[6] = mesh.add_vertex(MyMesh::Point( 1,  1, -1));
    sommets[7] = mesh.add_vertex(MyMesh::Point(-1,  1, -1));


    // on construit des faces à partir des sommets

    std::vector<MyMesh::VertexHandle> uneNouvelleFace;

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[0]);
    uneNouvelleFace.push_back(sommets[1]);
    uneNouvelleFace.push_back(sommets[2]);
    uneNouvelleFace.push_back(sommets[3]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[7]);
    uneNouvelleFace.push_back(sommets[6]);
    uneNouvelleFace.push_back(sommets[5]);
    uneNouvelleFace.push_back(sommets[4]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[1]);
    uneNouvelleFace.push_back(sommets[0]);
    uneNouvelleFace.push_back(sommets[4]);
    uneNouvelleFace.push_back(sommets[5]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[2]);
    uneNouvelleFace.push_back(sommets[1]);
    uneNouvelleFace.push_back(sommets[5]);
    uneNouvelleFace.push_back(sommets[6]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[3]);
    uneNouvelleFace.push_back(sommets[2]);
    uneNouvelleFace.push_back(sommets[6]);
    uneNouvelleFace.push_back(sommets[7]);
    mesh.add_face(uneNouvelleFace);

    uneNouvelleFace.clear();
    uneNouvelleFace.push_back(sommets[0]);
    uneNouvelleFace.push_back(sommets[3]);
    uneNouvelleFace.push_back(sommets[7]);
    uneNouvelleFace.push_back(sommets[4]);
    mesh.add_face(uneNouvelleFace);

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh,1);

}

/* **** fin de la partie boutons et IHM **** */



/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}


// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, int widget_nbr, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }

    if(widget_nbr==1)
        ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);
    if(widget_nbr==2)
        ui->displayWidget2->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);
    if(widget_nbr==3)
        ui->displayWidget3->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);


    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    if(widget_nbr==1)
        ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);
    if(widget_nbr==2)
        ui->displayWidget2->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);
    if(widget_nbr==3)
        ui->displayWidget3->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    if(widget_nbr==1)
        ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);
    if(widget_nbr==2)
        ui->displayWidget2->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);
    if(widget_nbr==3)
        ui->displayWidget3->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}
