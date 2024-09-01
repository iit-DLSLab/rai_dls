#include <Gui/RenderData.h>

#include <Kin/kin.h>
#include <Kin/frame.h>
#include <Gui/opengl.h>
#include <Algo/spline.h>

int main( void ){
  rai::Configuration C;
  rai::Frame *f = C.addFile(rai::raiPath("../rai-robotModels/panda/panda.g"));
  f->set_X()->addRelativeRotationDeg(90,0,0,1);

//  C.view(false, "test");

  RenderScene scene;
//  scene.drawShadows=false;
//  scene.drawTransparents=false;

  // meshes of C
  for(rai::Frame *f:C.frames){
    if(f->shape && f->shape->type()==rai::ST_mesh){
      scene.add().mesh(f->shape->mesh(), f->ensure_X());
    }
  }

  { // floor
    rai::Mesh m;
    m.setQuad();
    m.scale(5., 5., 0.);
    m.C = {1., .95, .9};
    scene.add().mesh(m, 0);
  }

  scene.addAxes(.2, "t(.2 .2 .2)");

  scene.addText("bla:dat 0.13098 t sec() []", 10.0, 20.0, 1.);

  if(true){ // ball
    rai::Mesh m;
    m.setSSBox(.4, .4, .4, .1, 2);
    //m.setSphere(0); m.scale(.2);
    m.C = {1., .5, .5, .5};
    scene.add().mesh(m, rai::Transformation("t(0 0 1.)"), .9);

    m.setSSBox(.2, .2, .2, .05, 2);
    m.C = {.5, 1., .5, .5};
    scene.add().mesh(m, rai::Transformation("t(0 -.3 1.)"), .9);

    m.setSSBox(.2, .2, .2, .05, 2);
    m.C = {.5, .5, 1.};
    scene.add().mesh(m, rai::Transformation("t(-.4 -.4 .4)"), .9);
  }

  scene.addLight({5.,5.,5.}, {0.,0.,1.});
  scene.addLight({-5.,0.,5.}, {0.,0.,1.});

  byteA img;
  read_ppm(img, "../opengl/box.ppm",false);
  add_alpha_channel(img, 120);
  scene.addQuad(img, 20, 20, 60, 60);

  OpenGL gl;
  gl.camera.setDefault();
#if 1
  gl.add(scene);
  gl.update(true);
#else
  gl.add(scene);
  gl.renderInBack(300, 500); //offscreen rendering tested
  write_png(gl.captureImage, "z.png");
  gl.update(true);
#endif

  return 0;
}
