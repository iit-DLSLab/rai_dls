/*  ------------------------------------------------------------------
    Copyright (c) 2011-2020 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#pragma once

#include "../Geo/geo.h"
#include "../Core/graph.h"
#include "../Geo/mesh.h"
#include "../Geo/signedDistanceFunctions.h"

/* TODO:
 * replace the types by more fundamental:
 *  shapes: ssbox or ssmesh -- nothing else
 *  joint: 7bits
 */

namespace rai {
struct Configuration;
struct Frame;
struct Dof;
struct Joint;
struct Shape;
struct Inertia;
struct ForceExchange;
struct ParticleDofs;
struct PathDof;
enum JointType { JT_none=0, JT_hingeX, JT_hingeY, JT_hingeZ, JT_transX, JT_transY, JT_transZ, JT_transXY, JT_trans3, JT_transXYPhi, JT_transYPhi, JT_universal, JT_rigid, JT_quatBall, JT_phiTransXY, JT_XBall, JT_free, JT_generic, JT_tau, JT_path };
enum BodyType  { BT_none=-1, BT_dynamic=0, BT_kinematic, BT_static, BT_soft };
}

typedef rai::Array<rai::Frame*> FrameL;
typedef rai::Array<rai::Joint*> JointL;
typedef rai::Array<rai::Dof*> DofL;
typedef rai::Array<rai::Shape*> ShapeL;

extern rai::Frame& NoFrame;
//extern rai::Shape& NoShape;
//extern rai::Joint& NoJoint;

namespace rai {

//===========================================================================
//
// Transformation tokens allow setting frame state (in destructor) after manipulating a transform
//

struct Transformation_Xtoken {
  Frame& f;
  Transformation_Xtoken(Frame& _f):f(_f) {}
  ~Transformation_Xtoken();
  Transformation* operator->();
  Transformation& operator*();
  void operator=(const Transformation&);
};

struct Transformation_Qtoken {
  Frame& f;
  Transformation_Qtoken(Frame& _f):f(_f) {}
  ~Transformation_Qtoken();
  Transformation* operator->();
  Transformation& operator*();
  void operator=(const Transformation&);
};

//===========================================================================

/// a Frame can have a link (also joint), shape (visual or coll), and/or intertia (mass) attached to it
struct Frame : NonCopyable {
  Configuration& C;        ///< a Frame is uniquely associated with a Configuration
  uint ID;                 ///< unique identifier (index in Configuration.frames)
  String name;             ///< name
  Frame* parent=nullptr;   ///< parent frame
  FrameL children;         ///< list of children
  Frame* prev=0;           ///< same frame in the previous time slice - if time sliced

 protected:
  Transformation Q=0;        ///< relative transform to parent
  Transformation X=0;        ///< frame's absolute pose
  //data structure state (lazy evaluation leave the state structure out of sync)
  bool _state_X_isGood=true; // X represents the current state
  void _state_setXBadinBranch();
  void _state_updateAfterTouchingX();
  void _state_updateAfterTouchingQ();
  //low-level fwd kinematics computation
  void calc_X_from_parent();
  void calc_Q_from_parent(bool enforceWithinJoint = true);

 public:
  double tau=0.;              ///< frame's relative time transformation (could be thought as part of the transformation X in space-time)
  std::shared_ptr<Graph> ats; ///< list of any-type attributes

  //attachments to the frame
  Joint* joint=nullptr;          ///< this frame is an articulated joint
  Shape* shape=nullptr;          ///< this frame has a (collision or visual) geometry
  Inertia* inertia=nullptr;      ///< this frame has inertia (is a mass)
  //TODO have a single list of all attached dofs (also joint)
  Array<ForceExchange*> forces;  ///< this frame exchanges forces with other frames
  ParticleDofs* particleDofs=nullptr; ///< this frame is a set of particles that are dofs themselves
  PathDof* pathDof=nullptr; ///< this frame is a set of particles that are dofs themselves

  Frame(Configuration& _C, const Frame* copyFrame=nullptr);
  Frame(Frame* _parent);
  ~Frame();

  //accessors to attachments
  Shape& getShape();
  Inertia& getInertia();

  //accessors to transformations
  const Transformation& ensure_X();
  const Transformation& get_Q() const;
  const Transformation& get_X() const;
  Transformation_Xtoken set_X() { return Transformation_Xtoken(*this); }
  Transformation_Qtoken set_Q();

  //structural operations
  Frame& setParent(Frame* _parent, bool keepAbsolutePose_and_adaptRelativePose=false, bool checkForLoop=false);
  Frame& unLink();
  Frame* insertPreLink(const rai::Transformation& A=0);
  Frame* insertPostLink(const rai::Transformation& B=0);

  //structural information/retrieval
  bool isChildOf(const Frame* par, int order=1) const;
  void getRigidSubFrames(FrameL& F, bool includeRigidJoints=false) const; ///< recursively collect all rigidly attached sub-frames (e.g., shapes of a link), (THIS is not included)
  void getPartSubFrames(FrameL& F) const; ///< recursively collect all frames of this part
  void getSubtree(FrameL& F) const;
  FrameL getSubtree() const { FrameL F; getSubtree(F); return F; }
  Frame* getRoot();
  FrameL getPathToRoot();
  Frame* getUpwardLink(rai::Transformation& Qtotal=NoTransformation, bool untilPartBreak=false) const; ///< recurse upward BEFORE the next joint and return relative transform (this->Q is not included!b)
  Frame* getDownwardLink(bool untilPartBreak=false) const; ///< recurse upward BEFORE the next joint and return relative transform (this->Q is not included!b)
  FrameL getPathToUpwardLink(bool untilPartBreak=false); ///< recurse upward BEFORE the next joint and return relative transform (this->Q is not included!b)
  const char* isPart() const;
  Dof* getDof() const;

  void prefixSubtree(const char* prefix);
  void transformToDiagInertia();

  //composed object manipulation
  void computeCompoundInertia(bool clearChildInertias=true);
  void convertDecomposedShapeToChildFrames();

  //I/O
  void read(const Graph& ats);
  void write(Graph& G);
  void write(std::ostream& os) const;

  //-- HIGHER LEVEL USER INTERFACE
  Frame& setShape(rai::ShapeType shape, const arr& size);
  Frame& setPose(const rai::Transformation& _X);
  Frame& setPosition(const arr& pos);
  Frame& setQuaternion(const arr& quat);
  Frame& setRelativePose(const rai::Transformation& _Q);
  Frame& setRelativePosition(const arr& pos);
  Frame& setRelativeQuaternion(const arr& quat);
  Frame& setPointCloud(const arr& points, const byteA& colors={});
  Frame& setConvexMesh(const arr& points, const byteA& colors={}, double radius=0.);
  Frame& setMesh(const rai::Mesh& m);
  Frame& setSdf(const SDF_GridData& sdf);
  Frame& setColor(const arr& color);
  Frame& setJoint(rai::JointType jointType);
  Frame& setContact(int cont);
  Frame& setMass(double mass);
  Frame& setAttribute(const char* key, double value);
  Frame& setJointState(const arr& q); ///< throws error if this frame is not also a joint, and if q.size() != joint->dim

  arr getPose() { return ensure_X().getArr7d(); }
  arr getPosition() { return ensure_X().pos.getArr(); }
  arr getQuaternion() { return ensure_X().rot.getArr4d(); }
  arr getRotationMatrix() { return ensure_X().rot.getArr(); }
  arr getRelativePosition() const { return get_Q().pos.getArr(); }
  arr getRelativeQuaternion() const { return get_Q().rot.getArr4d(); }
  arr getSize() const ;
  arr getMeshPoints() const ;
  uintA getMeshTriangles() const ;
  arr getMeshCorePoints() const ;
  arr getJointState() const; ///< throws error if this frame is not also a joint

  friend struct Configuration;
  friend struct Configuration_ext;
  friend struct KinematicSwitch;
  friend struct Joint;
  friend struct Transformation_Xtoken;
  friend struct Transformation_Qtoken;
  friend struct ConfigurationViewer;
};
stdOutPipe(Frame)

//===========================================================================

struct Dof {
  Frame* frame=0;    ///< this is the frame that this dof articulates! I.e., the output frame
  bool active=true;  ///< if false, this dof is not considered part of the configuration's q-vector
  uint dim=UINT_MAX;
  uint qIndex=UINT_MAX;
  arr  limits;       ///< joint limits (lo, up, [maxvel, maxeffort])
  Joint* mimic=0;    ///< if non-nullptr, this joint's state is identical to another's
  JointL mimicers;   ///< list of mimicing joints
  bool isStable=false;

  // sampling info:
  double sampleUniform=0.; //prob for uniform initialization within limits
  double sampleSdv=.01; //sdv of gaussian around default
  arr q0; //mean of gaussian, if not defined -> copyX from prev

  virtual ~Dof() {}
  virtual void setDofs(const arr& q, uint n=0) = 0;
  virtual arr calcDofsFromConfig() const = 0;
  virtual void setRandom(uint timeSlices_d1, int verbose);
  arr getDofState();
  virtual String name() const = 0;

  void copyParameters(Dof* copy){
    qIndex=copy->qIndex; dim=copy->dim;
    limits=copy->limits; q0=copy->q0;
    active=copy->active;
    isStable=copy->isStable;
    sampleUniform=copy->sampleUniform;  sampleSdv=copy->sampleSdv;
  }
  void setActive(bool _active);

  const Joint* joint() const;
  const ForceExchange* fex() const;

  virtual void write(std::ostream& os) const;
};
stdOutPipe(Dof)

//===========================================================================

/// for a Frame with Joint-Link, the relative transformation 'Q' is articulated
struct Joint : Dof, NonCopyable {
  // joint information
  //  byte generator;    ///< (7bits), h in Featherstone's code (indicates basis vectors of the Lie algebra, but including the middle quaternion w)
  String code;       ///< for JT_generic: code "txyzwabc" to indicate transformations; dim==code.N
  double H=1.;       ///< control cost scalar
  double scale=1.;   ///< scaling robot-q = scale * q-vector

  Vector axis=0;          ///< joint axis (same as X.rot.getX() for standard hinge joints)
  Enum<JointType> type;   ///< joint type

  //attachments to the joint
  //struct Uncertainty* uncertainty=nullptr;

  Joint(Frame& f, JointType type);
  Joint(Frame& f, Joint* copyJoint=nullptr);
  Joint(Frame& from, Frame& f, Joint* copyJoint=nullptr);
  ~Joint();

  const Transformation& X() const; ///< the base frame, where the joint STARTS (i.e. parent->X)
  const Transformation& Q() const; ///< the transformation realized by this joint (i.e. from parent->X to frame->X)
  Frame* from() const { return frame->parent; }
  virtual String name() const { return STRING(frame->name<<'.'<<frame->ID); }

  void setMimic(Joint* j, bool unsetPreviousMimic=false);
  void setDofs(const arr& q, uint n=0);
  arr calcDofsFromConfig() const;
  arr getScrewMatrix();
  uint getDimFromType() const;
  arr get_h() const;

  bool isPartBreak();

  //access the C's q vector
  double& get_q();

  //structural operations
  void makeRigid();
  void makeFree(double H_cost=0.);
  void setType(JointType _type);
  void setGeneric(const char* _code);
  void flip();

  void read(const Graph& ats);
  void write(Graph& g);
  void write(std::ostream& os) const;
};
stdOutPipe(Joint)

//===========================================================================

/// a Frame with Inertia has mass and, in physical simulation, has forces associated with it
struct Inertia : NonCopyable {
  Frame& frame;
  double mass=-1.;
  Matrix matrix=0;
  Enum<BodyType> type;
  Vector com=0;             ///< its center of mass

  Inertia(Frame& f, rai::Inertia* copyInertia=nullptr);
  ~Inertia();

  void setZero(){ mass=0; com=0; matrix=0; }
  void add(const Inertia& I, const rai::Transformation& rel);
  void defaultInertiaByShape();

  void write(std::ostream& os) const;
  void write(Graph& g);
  void read(const Graph& G);
};
stdOutPipe(Inertia)

//===========================================================================

/// a Frame with Shape is a collision or visual object
struct Shape : NonCopyable, GLDrawer {
  Frame& frame;
  Enum<ShapeType> _type;
  arr size;
  shared_ptr<Mesh> _mesh;
  shared_ptr<Mesh> _sscCore;
  shared_ptr<SDF_GridData> _sdf;
  char cont=0;           ///< are contacts registered (or filtered in the callback)

  double radius() { if(size.N) return size(-1); return 0.; }
  Enum<ShapeType>& type() { return _type; }
  Mesh& mesh() { if(!_mesh){ if(_type==ST_none) _type=ST_mesh; _mesh = make_shared<Mesh>(); } return *_mesh; }
  Mesh& sscCore() { if(!_sscCore){ if(_type==ST_none) _type=ST_ssCvx;  _sscCore = make_shared<Mesh>();  } return *_sscCore; }
  SDF_GridData& sdf() { if(!_sdf){ if(_type==ST_none) _type=ST_sdf; _sdf = make_shared<SDF_GridData>(); } return *_sdf; }
  double alpha() { arr& C=mesh().C; if(C.N==4 || C.N==2) return C(-1); return 1.; }

  void createMeshes();
  shared_ptr<ScalarFunction> functional(bool worldCoordinates=true);

  Shape(Frame& f, const Shape* copyShape=nullptr); //new Shape, being added to graph and frame's shape lists
  virtual ~Shape();

  bool canCollideWith(const Frame* f) const;

  void read(const Graph& ats);
  void write(std::ostream& os) const;
  void write(Graph& g);
  void glDraw(OpenGL&);
};

//===========================================================================

}// namespace rai

