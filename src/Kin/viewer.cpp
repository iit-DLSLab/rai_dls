/*  ------------------------------------------------------------------
    Copyright (c) 2011-2020 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#include "viewer.h"
#include "frame.h"
#include "../Gui/opengl.h"
#include <iomanip>

rai::ConfigurationViewer::~ConfigurationViewer() {
  if(gl) gl.reset();
}

OpenGL& rai::ConfigurationViewer::ensure_gl() {
  if(!gl) {
    gl = make_shared<OpenGL>("ConfigurationViewer");
    gl->camera.setDefault();
    gl->add(*this);
  }
  return *gl;
}

void rai::ConfigurationViewer::setCamera(rai::Frame* camF){
  ensure_gl();
  rai::Camera& cam = gl->camera;
  {
    auto _dataLock = gl->dataLock(RAI_HERE);
    if(camF){
      cam.X = camF->ensure_X();

      rai::Node *at=0;
      if((at=camF->ats->getNode("focalLength"))) cam.setFocalLength(at->as<double>());
      if((at=camF->ats->getNode("orthoAbsHeight"))) cam.setHeightAbs(at->as<double>());
      if((at=camF->ats->getNode("zRange"))){ arr z=at->as<arr>(); cam.setZRange(z(0), z(1)); }
      if((at=camF->ats->getNode("width"))) gl->width=at->as<double>();
      if((at=camF->ats->getNode("height"))) gl->height=at->as<double>();
      //    cam.setWHRatio((double)gl->width/gl->height);
    }else{
      gl->camera.setDefault();
    }
  }
  gl->resize(gl->width, gl->height);
}

int rai::ConfigurationViewer::update(const char* text, bool nonThreaded) { ensure_gl(); return gl->update(text, nonThreaded); }

int rai::ConfigurationViewer::watch(const char* text) { ensure_gl(); return gl->watch(text); }

void rai::ConfigurationViewer::add(GLDrawer& c) { ensure_gl(); gl->add(c); }

void rai::ConfigurationViewer::resetPressedKey() { ensure_gl(); gl->pressedkey=0; }

void rai::ConfigurationViewer::clear(){
  auto _dataLock = gl->dataLock(RAI_HERE);
  C.clear();
  framePath.clear();
}

int rai::ConfigurationViewer::update(bool watch) {
  ensure_gl();

  {
    auto _dataLock = gl->dataLock(RAI_HERE);
    gl->text = drawText;
    if(watch) gl->text <<"\n[press key]";
  }

  int ret=0;
  if(watch) {
    gl->raiseWindow();
    ret = gl->watch();
    gl->text = drawText;
  }else{
    ret = gl->update(nullptr, false);
  }

  return ret;
}

void rai::ConfigurationViewer::raiseWindow(){
  ensure_gl();
  gl->raiseWindow();
}

int rai::ConfigurationViewer::setConfiguration(const rai::Configuration& _C, const char* text, bool watch) {
  bool copyMeshes = false;
  if(_C.frames.N!=C.frames.N) copyMeshes = true;
  else{
    for(uint i=0;i<C.frames.N;i++){
      rai::Shape *s = _C.frames.elem(i)->shape;
      rai::Shape *r = C.frames.elem(i)->shape;
      if((!s) != (!r)){ copyMeshes=true; break; }
      if(!s) continue;
      if(s->_type != r->_type){ copyMeshes=true; break; }
      if(s->size != r->size){ copyMeshes=true; break; }
      if(s->_mesh && r->_mesh && (s->_mesh.get() != r->_mesh.get())){ copyMeshes=true; break; }
      if(s->_mesh && s->_mesh->listId<0){ copyMeshes=true; break; }
    }
  }
  if(copyMeshes) recopyMeshes(_C);

  CHECK_EQ(_C.frames.N, C.frames.N, "");

  ensure_gl();

  if(_C.proxies.N) {
    auto _dataLock = gl->dataLock(RAI_HERE);
    C.copyProxies(_C.proxies);
  }

  {
    auto _dataLock = gl->dataLock(RAI_HERE);
    for(rai::Frame* f:_C.frames) CHECK(f->_state_X_isGood, "");
#if 0
    framePath = _C.getFrameState();
#else
    framePath.resize(_C.frames.N, 7).setZero();
    for(uint i=0;i<framePath.d0;i++){
      rai::Frame *f = _C.frames.elem(i);
      if(f->shape) framePath[i] = f->ensure_X().getArr7d();
    }
#endif
    framePath.reshape(1, _C.frames.N, 7);
    drawTimeSlice=0;
    drawSubFrames.clear();
    if(text) drawText = text;
  }

  return update(watch);
}

int rai::ConfigurationViewer::setPath(ConfigurationL& Cs, const char* text, bool watch) {
  CHECK(C.frames.N, "setPath requires that you setConfiguration first");

  uintA frames;
  frames.setStraightPerm(Cs.first()->frames.N);

  arr X(Cs.N, frames.N, 7);
  for(uint t=0; t<X.d0; t++) {
    for(uint i=0; i<X.d1; i++) {
      X(t, i, {}) = Cs(t)->frames(frames(i))->getPose();
    }
  }

  return setPath(X, text, watch);
}

int rai::ConfigurationViewer::setPath(rai::Configuration& _C, const arr& jointPath, const char* text, bool watch, bool full) {
  setConfiguration(_C, 0, false);
  CHECK(C.frames.N, "setPath requires that you setConfiguration first");

  arr X(jointPath.d0, _C.frames.N, 7);
  for(uint t=0; t<X.d0; t++) {
    _C.setJointState(jointPath[t]);
    for(uint i=0; i<X.d1; i++) {
      X(t, i, {}) = _C.frames.elem(i)->getPose();
    }
  }

  return setPath(X, text, watch);
}

int rai::ConfigurationViewer::setPath(const arr& _framePath, const char* text, bool watch, bool full) {
  CHECK(C.frames.N, "setPath requires that you setConfiguration first");

  CHECK_EQ(_framePath.nd, 3, "");
  CHECK_EQ(_framePath.d1, C.frames.N, "");
  CHECK_EQ(_framePath.d2, 7, "");

  {
    auto _dataLock = gl->dataLock(RAI_HERE);
    framePath = _framePath;
    drawFullPath = full;
    drawTimeSlice=-1;
    if(text) drawText = text;
  }

  return update(watch);
}

bool rai::ConfigurationViewer::playVideo(uint T, uint nFrames, bool watch, double delay, const char* saveVideoPath){
  if(rai::getDisableGui()) return false;

  const rai::String tag = drawText;

  if(saveVideoPath) {
    rai::system(STRING("mkdir -p " <<saveVideoPath));
    rai::system(STRING("rm -f " <<saveVideoPath <<"*.png"));
  }

  CHECK_GE(C.frames.N, T*nFrames, "");

  FrameL F;
  if(C.frames.nd==1) F = C.frames({0,T*nFrames-1});
  else F = C.frames;
  F.reshape(T, nFrames);

  for(uint t=0; t<F.d0; t++) {
    {
      auto _dataLock = gl->dataLock(RAI_HERE);
      drawSubFrames = F[t];
      drawText.clear() <<tag <<"\n(slice " <<t <<'/' <<F.d0;
      if(phaseFactor>0.) drawText <<", phase " <<phaseFactor*(double(t)+phaseOffset);
      drawText <<")";
    }

    if(delay<0.) {
      update(true);
    } else {
      update(false);
      if(delay) rai::wait(delay / F.d0);
    }

    {
      auto _dataLock = gl->dataLock(RAI_HERE);
      if(saveVideoPath) write_png(gl->captureImage, STRING(saveVideoPath<<std::setw(4)<<std::setfill('0')<<t<<".png"));
    }
  }
  drawText = tag;
  if(watch && rai::getInteractivity()) {
    int key = update(true);
    drawSubFrames.clear();
    return !(key==27 || key=='q');
  }
  drawSubFrames.clear();
  return false;

}

bool rai::ConfigurationViewer::playVideo(bool watch, double delay, const char* saveVideoPath) {
#ifndef RAI_GL
  return false;
#endif
  const rai::String tag = drawText;

  if(saveVideoPath) {
    rai::system(STRING("mkdir -p " <<saveVideoPath));
    rai::system(STRING("rm -f " <<saveVideoPath <<"*.png"));
  }

  for(uint t=0; t<framePath.d0; t++) {
    {
      auto _dataLock = gl->dataLock(RAI_HERE);
      drawTimeSlice=t;
      drawText.clear() <<tag <<" (config:" <<t <<'/' <<framePath.d0 <<")";
    }

    if(delay<0.) {
      update(true);
    } else {
      update(false);
      if(delay) rai::wait(delay / framePath.d0);
    }

    {
      auto _dataLock = gl->dataLock(RAI_HERE);
      if(saveVideoPath) write_png(gl->captureImage, STRING(saveVideoPath<<std::setw(4)<<std::setfill('0')<<t<<".png"));
    }
  }
  drawText = tag;
  if(watch && rai::getInteractivity()) {
    int key = update(true);
    return !(key==27 || key=='q');
  }
  return false;
}

void rai::ConfigurationViewer::savePng(const char* saveVideoPath) {
  write_png(gl->captureImage, STRING(saveVideoPath<<std::setw(4)<<std::setfill('0')<<(pngCount++)<<".png"), true);
}

rai::Camera& rai::ConfigurationViewer::displayCamera() {
  ensure_gl();
  return gl->camera;
}

byteA rai::ConfigurationViewer::getRgb() {
  ensure_gl();
  byteA image = gl->captureImage;
  flip_image(image);
  return image;
}

floatA rai::ConfigurationViewer::getDepth() {
  ensure_gl();
  floatA depth = gl->captureDepth;
  flip_image(depth);
  for(float& d:depth) {
    if(d==1.f || d==0.f) d=-1.f;
    else d = gl->camera.glConvertToTrueDepth(d);
  }
  return depth;
}

void rai::ConfigurationViewer::recopyMeshes(const rai::Configuration& _C) {
  ensure_gl();

  {
    gl->dataLock.lock(RAI_HERE);
    C.copy(_C, false);
    gl->dataLock.unlock();
    //deep copy meshes!
//    for(rai::Frame* f:C.frames) if(f->shape) {
//        shared_ptr<Mesh> org = f->shape->_mesh;
//        if(org){
//          f->shape->_mesh = make_shared<Mesh> (*org.get());
//          f->shape->_mesh->listId=0;
//        }
//      }
  }
  rai::Frame *camF = C.getFrame("camera_gl", false);
  if(camF) setCamera(camF);
}

void rai::ConfigurationViewer::glDraw(OpenGL& gl) {
#ifdef RAI_GL
  glStandardScene(NULL, gl);

  if(!framePath.N) {
    gl.text <<"\nConfigurationViewer: NOTHING TO DRAW";
    return;
  }

  glPushMatrix();

  rai::Transformation T;

  //draw frame paths
  if(drawFrameLines) {
    glColor(0., 0., 0., .2);
    glLoadIdentity();
    for(uint i=0; i<framePath.d1; i++) {
      glBegin(GL_LINE_STRIP);
      for(uint t=0; t<framePath.d0; t++) {
        T.set(&framePath(t, i, 0));
        //          glTransform(pose);
        glVertex3d(T.pos.x, T.pos.y, T.pos.z);
      }
      glEnd();
    }
  }

  if(drawSubFrames.N){
//    C.setFrameState(framePath[t]);
    C.glDraw_sub(gl, drawSubFrames, 0);
  } else if(drawTimeSlice>=0) {
    uint t=drawTimeSlice;
    CHECK_LE(t+1, framePath.d0, "");
    CHECK_EQ(framePath.d1, C.frames.N, "");
    CHECK_EQ(framePath.d2, 7, "");

    C.setFrameState(framePath[t]);
    C.glDraw_sub(gl, C.frames, 0);
  } else {
    if(drawFullPath) {
      CHECK_EQ(framePath.d1, C.frames.N, "");
      CHECK_EQ(framePath.d2, 7, "");
      for(uint t=0; t<framePath.d0; t++) {
        C.setFrameState(framePath[t]);
        C.glDraw_sub(gl, C.frames, 1); //opaque
      }
      for(uint t=0; t<framePath.d0; t++) {
        C.setFrameState(framePath[t]);
        C.glDraw_sub(gl, C.frames, 2); //transparent
      }
    } else {
      NIY;
    }
  }

  glPopMatrix();
#else
  NICO
#endif
}

