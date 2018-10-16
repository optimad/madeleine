#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_CG.hpp"
#include "surfunstructured.hpp"
#include "commons.hpp"
#include <exception>

using namespace bitpit;

// =================================================================================== //

void test00001() {


    SurfUnstructured sphere1(2,3);
    sphere1.setExpert(true);
    sphere1.importSTL("./data/sphere1.stl", false);

    PiercedVector<Vertex> vv1;
    sphere1.deleteCoincidentVertices();
    vv1 = sphere1.getVertices();

    std::cout << " n vertices : " << sphere1.getVertexCount() << std::endl;
    std::cout << " n cells : " << sphere1.getInternalCount() << std::endl;

    PiercedVector<double> data1;
    for (auto v : vv1){
        data1.insert(v.getId(), std::sin(v[0]));
    }

    {
        vector<double> vdata1;
        for (auto v : vv1){
            vdata1.push_back(std::sin(v[0]));
        }
        sphere1.getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::POINT, vdata1) ;
        sphere1.getVTK().setName("sphere1");
        sphere1.write();
    }

    SurfaceSkdTree tree1(&sphere1);
    tree1.build();

    SurfUnstructured sphere2(2,3);
    sphere2.setExpert(true);
    sphere2.importSTL("./data/sphere2.stl", false);

    PiercedVector<Vertex> vv2;
    sphere2.deleteCoincidentVertices();
    vv2 = sphere2.getVertices();

    std::cout << " n vertices : " << sphere2.getVertexCount() << std::endl;
    std::cout << " n cells : " << sphere2.getInternalCount() << std::endl;

    vector<double> dataRef;

    for (auto v : vv2){
        dataRef.push_back(std::sin(v[0]));
    }

    sphere2.getVTK().setName("sphere2");
    sphere2.getVTK().addData("dataRef", VTKFieldType::SCALAR, VTKLocation::POINT, dataRef) ;

    PiercedVector<double> data2;
    for (auto v : vv2){
        darray3 x = v.getCoords();
        long id;
        double dist;
        tree1.findPointClosestCell(x, &id, &dist);

        ConstProxyVector<long> vIds = sphere1.getCell(id).getVertexIds();

        darray3 lambda;
        darray3 xP = CGElem::projectPointTriangle(x, vv1[vIds[0]].getCoords(), vv1[vIds[1]].getCoords(), vv1[vIds[2]].getCoords(), lambda);

        double data = 0.;
        data += data1[vIds[0]]*lambda[0];
        data += data1[vIds[1]]*lambda[1];
        data += data1[vIds[2]]*lambda[2];

        data2.insert(v.getId(), data);

    }



    vector<double> vdata2;
    for (auto v : vv2){
        vdata2.push_back(data2[v.getId()]);
    }

    sphere2.getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::POINT, vdata2) ;
    sphere2.write();



    return;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if ENABLE_MPI==1
    MPI::Init(argc, argv);

    {
#endif
        /**<Calling mimmo Test routines*/
        try{
            test00001() ;
        }
        catch(std::exception & e){
            std::cout<<"example_00001 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if ENABLE_MPI==1
    }

    MPI::Finalize();
#endif

    return 0;

}
