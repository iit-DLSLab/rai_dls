/*  ------------------------------------------------------------------
    Copyright (c) 2011-2020 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#include "F_qFeatures.h"
#include "frame.h"
#include <climits>

//===========================================================================

F_qItself::F_qItself(bool relative_q0) : relative_q0(relative_q0) {}

F_qItself::F_qItself(PickMode pickMode, const StringA &picks, const rai::Configuration &C, bool relative_q0)
    : relative_q0(relative_q0)
{
  if (pickMode == allActiveJoints)
  {
    for (rai::Frame *f : C.frames)
      if (f->joint && f->parent && f->joint->active && f->joint->dim != 0)
      {
        frameIDs.append(f->ID);
        frameIDs.append(f->parent->ID);
      }
    frameIDs.reshape(-1, 2);
  }
  else if (pickMode == byJointNames)
  {
    for (rai::String s : picks)
    {
      if (s(-2) == ':')
        s.resize(s.N - 2, true);
      rai::Frame *f = C.getFrame(s);
      if (!f)
        HALT("pick '" << s << "' not found");
      if (!f->joint)
        HALT("pick '" << s << "' is not a joint");
      frameIDs.setAppend(f->ID);
    }
  }
  else if (pickMode == byExcludeJointNames)
  {
    for (rai::Joint *j : C.activeJoints)
    {
      if (picks.contains(j->frame->name))
        continue;
      frameIDs.setAppend(j->frame->ID);
    }
  }
  else
  {
    NIY
  }
}

F_qItself::F_qItself(const uintA &_selectedFrames, bool relative_q0)
    : relative_q0(relative_q0)
{
  frameIDs = _selectedFrames;
  fs = FS_qItself;
}

void F_qItself::phi(arr &q, arr &J, const rai::Configuration &C)
{
  CHECK(C._state_q_isGood, "");
  if (!frameIDs.nd)
  {
    q = C.getJointState();
    if (relative_q0)
    {
      for (rai::Joint *j : C.activeJoints)
        if (j->q0.N && j->qDim() == 1)
          q(j->qIndex) -= j->q0.scalar();
    }
    if (!!J)
      J.setId(q.N);
  }
  else
  {
    uint n = dim_phi(C);
    C.kinematicsZero(q, J, n);
    uint m = 0;
    if (frameIDs.nd)
    {
      for (uint i = 0; i < frameIDs.d0; i++)
      {
        rai::Joint *j = 0;
        bool flipSign = false;
        if (frameIDs.nd == 1)
        {
          rai::Frame *f = C.frames.elem(frameIDs.elem(i));
          j = f->joint;
          CHECK(j, "selected frame " << frameIDs.elem(i) << " ('" << f->name << "') is not a joint");
        }
        else
        {
          rai::Frame *a = C.frames.elem(frameIDs(i, 0));
          rai::Frame *b = C.frames.elem(frameIDs(i, 1));
          if (a->parent == b)
            j = a->joint;
          else if (b->parent == a)
          {
            j = b->joint;
            flipSign = true;
          }
          else
            HALT("a and b are not linked");
          CHECK(j, "");
        }
        for (uint k = 0; k < j->dim; k++)
        {
          if (j->active)
          {
            q.elem(m) = C.q.elem(j->qIndex + k);
          }
          else
          {
            q.elem(m) = C.qInactive.elem(j->qIndex + k);
          }
          if (flipSign)
            q.elem(m) *= -1.;
          if (relative_q0 && j->q0.N)
            q.elem(m) -= j->q0(k);
          if (!!J && j->active)
          {
            if (flipSign)
              J.elem(m, j->qIndex + k) = -1.;
            else
              J.elem(m, j->qIndex + k) = 1.;
          }
          m++;
        }
      }
      CHECK_EQ(n, m, "");
    }
  }

  std::cout << "Jjoint " << J << "\n";
}

void F_qItself::phi2(arr &q, arr &J, const FrameL &F)
{
  if (order != 0)
  {
    phi_finiteDifferenceReduce(q, J, F);
    return;
  }
  uint n = dim_phi2(F);
  if (!n)
  {
    q.clear();
    J.clear();
    return;
  }
  rai::Configuration &C = F.last()->C;
  CHECK(C._state_q_isGood, "");
  C.kinematicsZero(q, J, n);
  uint m = 0;
  CHECK(F.d0 == 1, "");
  FrameL FF = F[0];
  //  FF.reshape(-1,2);
  for (uint i = 0; i < FF.d0; i++)
  {
    rai::Joint *j = 0;
    bool flipSign = false;
    if (FF.nd == 1)
    {
      rai::Frame *f = FF.elem(i);
      j = f->joint;
      CHECK(j, "selected frame " << FF.elem(i) << " ('" << f->name << "') is not a joint");
    }
    else
    {
      rai::Frame *a = FF(i, 0);
      rai::Frame *b = FF(i, 1);
      if (a->parent == b)
        j = a->joint;
      else if (b->parent == a)
      {
        j = b->joint;
        flipSign = true;
      }
      else
        HALT("a and b are not linked");
      CHECK(j, "");
    }
    for (uint k = 0; k < j->dim; k++)
    {
      if (j->active)
      {
        // std::cout << "j->qIndex " << j->qIndex << "\n";
        // std::cout << "j->dim " << j->dim << "\n";

        q.elem(m) = C.q.elem(j->qIndex + k);
      }
      else
      {
        q.elem(m) = C.qInactive.elem(j->qIndex + k);
      }
      if (flipSign)
        q.elem(m) *= -1.;
      if (relative_q0 && j->q0.N)
        q.elem(m) -= j->q0(k);
      if (!!J && j->active)
      {
        if (flipSign)
          J.elem(m, j->qIndex + k) = -1.;
        else
          J.elem(m, j->qIndex + k) = 1.;
      }
      m++;
    }
  }
  CHECK_EQ(n, m, "");
}

uint F_qItself::dim_phi(const rai::Configuration &C)
{
  if (frameIDs.nd)
  {
    uint n = 0;
    for (uint i = 0; i < frameIDs.d0; i++)
    {
      rai::Joint *j = 0;
      if (frameIDs.nd == 1)
      {
        rai::Frame *f = C.frames.elem(frameIDs.elem(i));
        j = f->joint;
        CHECK(j, "selected frame " << frameIDs.elem(i) << " ('" << f->name << "') is not a joint");
      }
      else
      {
        rai::Frame *a = C.frames.elem(frameIDs(i, 0));
        rai::Frame *b = C.frames.elem(frameIDs(i, 1));
        if (a->parent == b)
          j = a->joint;
        else if (b->parent == a)
          j = b->joint;
        else
          HALT("a (" << a->name << ") and b (" << b->name << ") are not linked");
        CHECK(j, "");
      }
      n += j->qDim();
    }
    return n;
  }
  return C.getJointStateDimension();
}

uint F_qItself::dim_phi2(const FrameL &F)
{
  uint m = 0;
  FrameL FF = F[0];
  for (uint i = 0; i < FF.d0; i++)
  {
    rai::Joint *j = 0;
    if (FF.nd == 1)
    {
      rai::Frame *f = FF.elem(i);
      std::cout << "f name" << f->name << "\n";
      CHECK(j, "selected frame " << FF.elem(i) << " ('" << f->name << "') is not a joint");
    }
    else
    {
      rai::Frame *a = FF(i, 0);
      rai::Frame *b = FF(i, 1);
      if (a->parent == b)
        j = a->joint;
      else if (b->parent == a)
        j = b->joint;
      CHECK(j, "a (" << a->name << ") and b (" << b->name << ") are not linked");
    }
    m += j->dim;
  }
  return m;
}

//===========================================================================

void F_qZeroVel::phi2(arr &y, arr &J, const FrameL &F)
{
  CHECK_EQ(order, 1, "");
  F_qItself()
      .setOrder(order)
      .eval(y, J, F);
#if 1
  rai::Frame *f = F.last();
  if (f->joint->type == rai::JT_transXYPhi)
  {
    arr s = ARR(10., 10., 1.);
    y = s % y;
    if (!!J)
      J = s % J;
  }
  if (f->joint->type == rai::JT_free)
  {
    arr s = ARR(10., 10., 10., 1., 1., 1., 1.);
    y = s % y;
    if (!!J)
      J = s % J;
  }
#endif
}

uint F_qZeroVel::dim_phi2(const FrameL &F)
{
  return F_qItself()
      .setOrder(order)
      .dim(F);
}

//===========================================================================

rai::Array<rai::Joint *> getMatchingJoints(const ConfigurationL &Ktuple, bool zeroVelJointsOnly)
{
  rai::Array<rai::Joint *> matchingJoints;
  rai::Array<rai::Joint *> matches(Ktuple.N);
  bool matchIsGood;

  rai::Joint *j;
  for (rai::Frame *f : Ktuple.last()->frames)
    if ((j = f->joint) && j->active && !zeroVelJointsOnly)
    {
      matches.setZero();
      matches.last() = j;
      matchIsGood = true;

      for (uint k = 0; k < Ktuple.N - 1; k++)
      { // go through other configs
        if (Ktuple(k)->frames.N <= j->frame->ID)
        {
          matchIsGood = false;
          break;
        }
        rai::Frame *fmatch = Ktuple(k)->frames.elem(j->frame->ID);
        if (!fmatch)
        {
          matchIsGood = false;
          break;
        }
        rai::Joint *jmatch = fmatch->joint; // getJointByBodyIndices(j->from()->ID, j->frame->ID);
        if (!jmatch || j->type != jmatch->type)
        {
          matchIsGood = false;
          break;
        }
        if (j->from() && j->from()->ID != jmatch->from()->ID)
        {
          matchIsGood = false;
          break;
        }
        matches(k) = jmatch;
      }

      if (matchIsGood)
        matchingJoints.append(matches);
    }
  matchingJoints.reshape(matchingJoints.N / Ktuple.N, Ktuple.N);
  return matchingJoints;
}

//===========================================================================

void F_qLimits2::phi2(arr &y, arr &J, const FrameL &F)
{
  std::cout << "setting joint limits inside function " << "\n";

  uint M = dim_phi2(F);
  std::cout << "M " << M << "\n";
  std::cout << "LIMITS: J dimensions: " << J.d0 << " x " << J.d1 << "\n";

  F.last()->C.kinematicsZero(y, J, M);
  uint m = 0;
  for (rai::Frame *f : F)
  {
    rai::Joint *j = f->joint;
    if (!j)
      continue;
    if (!j->limits.N)
      continue;
    uint d = j->qDim();
    for (uint k = 0; k < d; k++)
    { // in case joint has multiple dimensions
      double lo = j->limits(2 * k + 0);
      double up = j->limits(2 * k + 1);
      uint i = j->qIndex + k;
      y.elem(m) = lo - f->C.q(i);
      if (!!J)
        J.elem(m, i) -= 1.;
      m++;
      y.elem(m) = f->C.q(i) - up;
      if (!!J)
        J.elem(m, i) += 1.;
      m++;
    }
  }
}

uint F_qLimits2::dim_phi2(const FrameL &F)
{
  uint m = 0;
  for (rai::Frame *f : F)
  {
    rai::Joint *j = f->joint;
    if (!j)
      continue;
    if (!j->limits.N)
      continue;
    m += 2 * j->qDim();
  }
  return m;
}

DofL getDofs(const FrameL &F)
{
  DofL dofs;
  for (rai::Frame *f : F)
  {
    if (f->joint && f->joint->active)
    {
      if (f->joint->limits.N)
        dofs.append(f->joint);
    }
    // for(rai::ForceExchange* fex:f->forces) if(&fex->a==f){
    //   if(fex->active && fex->limits.N) dofs.append(fex);
    // }
  }
  return dofs;
}

//===========================================================================
void F_qLimits::phi2(arr &y, arr &J, const FrameL &F)
{
  std::cout << "setting joint limits inside function " << "\n";
  uint M = dim_phi2(F);

  // std::cout << "M " << M << "\n";

  F.last()->C.kinematicsZero(y, J, M);
  CHECK(F.last()->C._state_q_isGood, "");
  uint m = 0;
  DofL dofs = getDofs(F);
  std::cout << "LIMITS: J dimensions: " << J.d0 << " x " << J.d1 << "\n";

  for (rai::Dof *dof : dofs)
    if (dof->limits.N)
    {
      // std::cout << "dof name " << *dof->name() << "\n";
      for (uint k = 0; k < dof->dim; k++)
      { // in case joint has multiple dimensions

        // std::cout << "dof dim " << dof->dim << "\n";
        // std::cout << "dof k " << k << "\n";

        double lo = dof->limits(2 * k + 0);
        double up = dof->limits(2 * k + 1);
        if (up >= lo)
        {
          uint i = dof->qIndex + k;
          // std::cout << "dof i " << i << "\n";

          double qi = F.last()->C.q(i);
          //        if(true){
          //          if(qi < lo) LOG(0) <<dof->name() <<' ' <<k <<' ' <<qi <<'<' <<lo <<" violates lower limit";
          //          if(qi > up) LOG(0) <<dof->name() <<' ' <<k <<' ' <<qi <<'>' <<up <<" violates upper limit";
          //        }
          y.elem(m) = lo - qi;
          if (!!J)
            J.elem(m, i) -= 1.;
          m++;
          y.elem(m) = qi - up;

          if (!!J)
            J.elem(m, i) += 1.;
          m++;
        }
        else
        {
          m += 2;
        }
      }
    }
  std::cout << "JLIMITS " << J << "\n";

  CHECK_EQ(m, M, "");
  // std::cout << "J " << J << "\n";
}

uint F_qLimits::dim_phi2(const FrameL &F)
{
  uint m = 0;
  DofL dofs = getDofs(F);
  for (rai::Dof *dof : dofs)
    if (dof->limits.N)
      m += 2 * dof->dim;
  // std::cout << "LIMITS JOINTS " << m << "\n";
  return m;
}

//===========================================================================

void F_qQuaternionNorms::phi2(arr &y, arr &J, const FrameL &F)
{
  uint n = dim_phi2(F);
  if (!n)
  {
    y.clear();
    J.clear();
    return;
  }
  rai::Configuration &C = F.first()->C;
  C.kinematicsZero(y, J, n);
  uint i = 0;
  for (const rai::Frame *f : F)
  {
    rai::Joint *j = f->joint;
    if (!j || !j->active)
      continue;
    if (j->type == rai::JT_quatBall || j->type == rai::JT_free || j->type == rai::JT_XBall)
    {
      arr q;
      if (j->type == rai::JT_quatBall)
        q.referToRange(C.q, j->qIndex + 0, j->qIndex + 3);
      if (j->type == rai::JT_XBall)
        q.referToRange(C.q, j->qIndex + 1, j->qIndex + 4);
      if (j->type == rai::JT_free)
        q.referToRange(C.q, j->qIndex + 3, j->qIndex + 6);
      double norm = sumOfSqr(q);
      y(i) = norm - 1.;

      if (!!J)
      {
        if (j->type == rai::JT_quatBall)
          for (uint k = 0; k < 4; k++)
            J.elem(i, j->qIndex + 0 + k) = 2. * q.elem(k);
        if (j->type == rai::JT_XBall)
          for (uint k = 0; k < 4; k++)
            J.elem(i, j->qIndex + 1 + k) = 2. * q.elem(k);
        if (j->type == rai::JT_free)
          for (uint k = 0; k < 4; k++)
            J.elem(i, j->qIndex + 3 + k) = 2. * q.elem(k);
      }
      i++;
    }
  }
}

uint F_qQuaternionNorms::dim_phi2(const FrameL &F)
{
  uint n = 0;
  for (const rai::Frame *f : F)
  {
    rai::Joint *j = f->joint;
    if (!j || !j->active)
      continue;
    if (j->type == rai::JT_quatBall || j->type == rai::JT_free || j->type == rai::JT_XBall)
      n++;
  }
  return n;
}

void F_qQuaternionNorms::setAllActiveQuats(const rai::Configuration &C)
{
  frameIDs.clear();
  for (const rai::Joint *j : C.activeJoints)
  {
    if (j->type == rai::JT_quatBall || j->type == rai::JT_free || j->type == rai::JT_XBall)
      frameIDs.append(j->frame->ID);
  }
}

//===========================================================================

void F_qTime::phi2(arr &y, arr &J, const FrameL &F)
{
  if (order == 0)
  {
    rai::Frame *f = F.scalar();
    double tau;
    f->C.kinematicsTau(tau, J, f);
    y.resize(1) = tau;
  }
  if (order == 1)
  { // WARNING: this is neg velocity... for ineq constraint
    CHECK_EQ(F.N, 2, "");
    arr y0, y1, J0, J1;
    order = 0;
    phi2(y0, J0, {F.elem(0)});
    phi2(y1, J1, {F.elem(1)});
    order = 1;
    y = y0 - y1;
    if (!!J)
      J = J0 - J1;
  }
  if (order == 2)
  {
    CHECK_EQ(F.N, 3, "");
    arr y0, y1, y2, J0, J1, J2;
    order = 0;
    phi2(y0, J0, {F.elem(0)});
    phi2(y1, J1, {F.elem(1)});
    phi2(y2, J2, {F.elem(2)});
    order = 2;
    y = y2 - 2. * y1 + y0;
    if (!!J)
      J = J2 - 2. * J1 + J0;
  }
}

//===========================================================================

// void F_feasibilityMargin::phi2(arr &y, arr &J, const FrameL &F)
// {

//   margin.states.linAcc = Eigen::Vector3d(0., 0., 0.);
//   margin.states.angAcc = Eigen::Vector3d(0., 0., 0.);
//   margin.states.extForce = Eigen::Vector3d(0., 0., 0.);
//   margin.states.extTorque = Eigen::Vector3d(0., 0., 0.);
//   // Convert to the base frame
//   margin.states.feetPos.resize(4);

//   // Set friction
//   margin.states.friction = 0.5;
//   // Set feet normals
//   margin.states.normals.resize(4);
//   margin.states.normals[0] = Eigen::Vector3d(0., 0., 1.);
//   margin.states.normals[1] = Eigen::Vector3d(0., 0., 1.);
//   margin.states.normals[2] = Eigen::Vector3d(0., 0., 1.);
//   margin.states.normals[3] = Eigen::Vector3d(0., 0., 1.);

//   uint M = dim_phi2(F);

//   F.last()->C.kinematicsZero(y, J, M);
//   std::cout << "MARGIN: J dimensions: " << J.d0 << " x " << J.d1 << "\n";

//   CHECK(F.last()->C._state_q_isGood, "");
//   uint m = 0;
//   DofL dofs = getDofs(F);
//   // for (int t = 0; t < 15; t++)
//   // {
//   for (rai::Dof *dof : dofs)
//   {
//     // std::cout << "dof is " << dof->name() << "\n";
//     // std::cout << "dof->qIndex is " << dof->qIndex << "\n";

//     // if ((dof->qIndex - 4) % 12 == 0 && dof->qIndex > 3)
//     if ((dof->qIndex - 4) % 12 == 0 && dof->qIndex > 12)
//     {

//       uint i = dof->qIndex;
//       // std::cout << "m index is" << m << "\n";
//       std::cout << "iIndex is" << i << "\n";

//       double pitch = F.last()->C.q(i);
//       double z = F.last()->C.q(i - 1);
//       double x = F.last()->C.q(i - 4);
//       double y_b = F.last()->C.q(i - 3);
//       double psi = F.last()->C.q(i - 2);

//       y.elem(m) = 0.15 - computeFeasibilityMargin(x, y_b, z, pitch, psi);

//       Transform inputTransform;

//       /// Prepare input in right frames
//       // Compute rotation matrix (from base to world) from roll/pitch/yaw - Intrinsic rotations.
//       Eigen::Matrix3d w_R_b = inputTransform.eulerToRotationMatrix(0.0, pitch, psi); // Make function of base rotation, pitch and yaw
//       Eigen::Matrix3d b_R_w = w_R_b.transpose();
//       margin.states.gravityRot = w_R_b.row(2);

//       // Convert to the base frame
//       margin.states.feetPos.resize(4);
//       // margin.states.feetPos[0] = b_R_w * (Eigen::Vector3d(0.44, 0.34, 0.0) - Eigen::Vector3d(x, y_b, z));   // LF Foot Position
//       // margin.states.feetPos[1] = b_R_w * (Eigen::Vector3d(0.44, -0.34, 0.0) - Eigen::Vector3d(x, y_b, z));  // RF Foot Position
//       // margin.states.feetPos[2] = b_R_w * (Eigen::Vector3d(-0.44, 0.34, 0.0) - Eigen::Vector3d(x, y_b, z));  // LH Foot Position
//       // margin.states.feetPos[3] = b_R_w * (Eigen::Vector3d(-0.44, -0.34, 0.0) - Eigen::Vector3d(x, y_b, z)); // RH Foot Position.
  
//       // margin.states.feetPos[0] = b_R_w * (Eigen::Vector3d(0.55, 0.34, 0.0) - Eigen::Vector3d(x, y_b, z));   // LF Foot Position
//       // margin.states.feetPos[1] = b_R_w * (Eigen::Vector3d(0.54, -0.32, 0.0) - Eigen::Vector3d(x, y_b, z));  // RF Foot Position
//       // margin.states.feetPos[2] = b_R_w * (Eigen::Vector3d(-0.33, 0.30, 0.0) - Eigen::Vector3d(x, y_b, z));  // LH Foot Position
//       // margin.states.feetPos[3] = b_R_w * (Eigen::Vector3d(-0.41, -0.4, 0.0) - Eigen::Vector3d(x, y_b, z)); // RH Foot Position.


//       margin.states.feetPos[0] = b_R_w * (Eigen::Vector3d(0.367, 0.32, 0.0) - Eigen::Vector3d(x, y_b, z));   // LF Foot Position
//       margin.states.feetPos[1] = b_R_w * (Eigen::Vector3d(0.47, -0.41, 0.0) - Eigen::Vector3d(x, y_b, z));  // RF Foot Position
//       margin.states.feetPos[2] = b_R_w * (Eigen::Vector3d(-0.44, 0.37, 0.0) - Eigen::Vector3d(x, y_b, z));  // LH Foot Position
//       margin.states.feetPos[3] = b_R_w * (Eigen::Vector3d(-0.41, -0.36, 0.0) - Eigen::Vector3d(x, y_b, z)); // RH Foot Position.


//       Eigen::MatrixXd footGradient = margin.getFeasibilityGradients(margin.states);
//       Eigen::Vector3d dMdG(footGradient(0), footGradient(1), footGradient(2));
//       Eigen::Vector3d dGdTheta(-cos(psi) * cos(pitch), cos(pitch) * sin(psi), -sin(pitch));
//       Eigen::Vector3d dGdPsi(sin(psi) * sin(pitch), sin(pitch) * cos(psi), 0.0);



//       J.elem(m, i - 4) = (footGradient(3) + footGradient(6) + footGradient(9) + footGradient(12)) * b_R_w(0, 0);
//       //  m++;

//       J.elem(m, i - 3) = (footGradient(4) + footGradient(7) + footGradient(10) + footGradient(13)) * b_R_w(1, 1);
//       // m++;

//       J.elem(m, i - 2) = -dMdG.dot(dGdPsi);
//       // m++;

//       J.elem(m, i - 1) = (footGradient(5) + footGradient(8) + footGradient(11) + footGradient(14)) * b_R_w(2, 2);
//       // m++;


//       std::cout << " margin.states.feetPos[0]  " <<  margin.states.feetPos[0]  << "\n";
//       std::cout << " margin.states.feetPos[1]  " <<  margin.states.feetPos[1] << "\n";
//       std::cout << " margin.states.feetPos[2]  " <<  margin.states.feetPos[2] << "\n";
//       std::cout << " margin.states.feetPos[3]  " <<  margin.states.feetPos[3] << "\n";

//       std::cout << "footGradient " << footGradient.transpose() << "\n";
//       std::cout << "b_R_w(2, 2) " << b_R_w(2, 2) << "\n";

//       std::cout << "Grad " << (footGradient(5) + footGradient(8) + footGradient(11) + footGradient(14)) * b_R_w(2, 2) << "\n";
//       std::cout << "z " << z << "\n";
//       std::cout << "x " << x << "\n";
//       std::cout << "y_b " << y_b << "\n";

//       std::cout << "margin " << computeFeasibilityMargin(x, y_b, z, pitch, psi) << "\n";

//       J.elem(m, i) = -dMdG.dot(dGdTheta);
//       double delta = 0.0175;
//       // auto marginPlusPitch = computeFeasibilityMargin(x, y_b, z, pitch + delta, psi);
//       // auto marginMinusPitch = computeFeasibilityMargin(x, y_b, z, pitch - delta, psi);
//       // J.elem(m, i) = 1000000 * (marginPlusPitch - marginMinusPitch);
//       // J.elem(m, i) = 1000000 * (marginPlusPitch - marginMinusPitch);

//       // m++;

//       // std::cout << "\n";

//       // std::cout << "gradPitch " << dMdG.dot(dGdTheta) << "\n";
//       // std::cout << "J.elem(m, i - 1) " << J.elem(m, i - 1) << "\n";

//       // std::cout << "gradPitch " << (marginPlusPitch - marginMinusPitch) << "\n";

//       // std::cout << "gradPitch " << J.elem(m, i) << "\n";
//       // std::cout << "pitch " << pitch << "\n";
//       // std::cout << "margin " << computeFeasibilityMargin(x, y_b, z, pitch, psi) << "\n";
//       // std::cout << "marginPlusPitch " << marginPlusPitch << "for pitch " << pitch + delta << "\n";
//       // std::cout << "marginNegatPitch " << marginMinusPitch << "for pitch " << pitch - delta << "\n";
//       // std::cout << "\n";

//       // Debug prints
//       // std::cout << "m: " << m << ", i: " << i << "\n";
//       // std::cout << "Assigned J(" << m << ", " << i - 4 << "), J(" << m << ", " << i - 3 << "), etc.\n";

//       // }
//       // }
//     }
//   }
//   std::cout << "JMargin " << J << "\n";
// }

// double F_feasibilityMargin::computeFeasibilityMargin(double x, double y, double z, double pitch, double psi)
// {
//   Transform inputTransform;

//   /// Prepare input in right frames
//   // Compute rotation matrix (from base to world) from roll/pitch/yaw - Intrinsic rotations.
//   Eigen::Matrix3d w_R_b = inputTransform.eulerToRotationMatrix(0.0, pitch, psi); // Make function of base rotation, pitch and yaw
//   Eigen::Matrix3d b_R_w = w_R_b.transpose();
//   margin.states.gravityRot = w_R_b.row(2);

//   // Convert to the base frame
//   margin.states.feetPos.resize(4);


//   // margin.states.feetPos[0] = b_R_w * (Eigen::Vector3d(0.55, 0.34, 0.0) - Eigen::Vector3d(x, y, z));   // LF Foot Position
//   // margin.states.feetPos[1] = b_R_w * (Eigen::Vector3d(0.54, -0.32, 0.0) - Eigen::Vector3d(x, y, z));  // RF Foot Position
//   // margin.states.feetPos[2] = b_R_w * (Eigen::Vector3d(-0.33, 0.30, 0.0) - Eigen::Vector3d(x, y, z));  // LH Foot Position
//   // margin.states.feetPos[3] = b_R_w * (Eigen::Vector3d(-0.41, -0.4, 0.0) - Eigen::Vector3d(x, y, z)); // RH Foot Position.

//       margin.states.feetPos[0] = b_R_w * (Eigen::Vector3d(0.367, 0.32, 0.0) - Eigen::Vector3d(x, y, z));   // LF Foot Position
//       margin.states.feetPos[1] = b_R_w * (Eigen::Vector3d(0.47, -0.41, 0.0) - Eigen::Vector3d(x, y, z));  // RF Foot Position
//       margin.states.feetPos[2] = b_R_w * (Eigen::Vector3d(-0.44, 0.37, 0.0) - Eigen::Vector3d(x, y, z));  // LH Foot Position
//       margin.states.feetPos[3] = b_R_w * (Eigen::Vector3d(-0.41, -0.36, 0.0) - Eigen::Vector3d(x, y, z)); // RH Foot Position.


//   // margin.states.feetPos[0] = b_R_w * (Eigen::Vector3d(0.38, 0.34, 0.2) - Eigen::Vector3d(x, y, z));   // LF Foot Position
//   // margin.states.feetPos[1] = b_R_w * (Eigen::Vector3d(0.46, -0.42, 0.2) - Eigen::Vector3d(x, y, z));  // RF Foot Position
//   // margin.states.feetPos[2] = b_R_w * (Eigen::Vector3d(-0.44, 0.37, 0.0) - Eigen::Vector3d(x, y, z));  // LH Foot Position
//   // margin.states.feetPos[3] = b_R_w * (Eigen::Vector3d(-0.41, -0.36, 0.0) - Eigen::Vector3d(x, y, z)); // RH Foot Position.

//   return margin.getFeasibilityMargin(margin.states);
// }

// void F_feasibilityMargin::init()
// {
//   std::string NewtorkPath = "/home/simple/ros_ws_lgp/src/reachability-aware-LGP/tamp_lgp/../../networks_minimal";
//   std::string parametersDirectory(NewtorkPath + "/examples/parameters/");

//   margin.loadDataScalersFromFile(parametersDirectory + "hyqreal/dataScalers.txt");
// }

// // Feasibility margin constraint

// void F_feasibilityMarginConstraint::phi2(arr &y, arr &J, const FrameL &F)
// {

//   margin.states.linAcc = Eigen::Vector3d(0., 0., 0.);
//   margin.states.angAcc = Eigen::Vector3d(0., 0., 0.);
//   margin.states.extForce = Eigen::Vector3d(0., 0., 0.);
//   margin.states.extTorque = Eigen::Vector3d(0., 0., 0.);
//   // Convert to the base frame
//   margin.states.feetPos.resize(4);

//   // Set friction
//   margin.states.friction = 0.5;
//   // Set feet normals
//   margin.states.normals.resize(4);
//   margin.states.normals[0] = Eigen::Vector3d(0., 0., 1.);
//   margin.states.normals[1] = Eigen::Vector3d(0., 0., 1.);
//   margin.states.normals[2] = Eigen::Vector3d(0., 0., 1.);
//   margin.states.normals[3] = Eigen::Vector3d(0., 0., 1.);

//   uint M = dim_phi2(F);

//   F.last()->C.kinematicsZero(y, J, M);
//   CHECK(F.last()->C._state_q_isGood, "");
//   uint m = 0;
//   DofL dofs = getDofs(F);
//   for (rai::Dof *dof : dofs)
//   {
//     // std::cout << "dof is " << dof << "\n";

//     if ((dof->qIndex - 4) % 12 == 0 && dof->qIndex > 3)
//     {
//       uint i = dof->qIndex;

//       double pitch = F.last()->C.q(i);
//       double z = F.last()->C.q(i - 1);
//       double x = F.last()->C.q(i - 4);
//       double y_b = F.last()->C.q(i - 3);
//       double psi = F.last()->C.q(i - 2);

//       double delta = 0.0175;
//       // auto marginPlusZ = computeFeasibilityMarginConstraint(x, y_b, z + delta, pitch, psi);
//       // auto marginMinusZ = computeFeasibilityMarginConstraint(x, y_b, z - delta, pitch, psi);

//       // auto marginPlusPitch = computeFeasibilityMarginConstraint(x, y_b, z, pitch + delta, psi);
//       // auto marginMinusPitch = computeFeasibilityMarginConstraint(x, y_b, z, pitch - delta, psi);

//       y.elem(m) = 0.05 - computeFeasibilityMarginConstraint(x, y_b, z, pitch, psi);

//       // y.elem(m) = 0.0;
//       // std::cout << "FM: " << computeFeasibilityMargin(x, y_b, z, pitch, psi) << "\n";
//       // std::cout << "CM: " << -100 * computeFeasibilityMargin(x, y_b, z, pitch, psi) << "\n";

//       // y.elem(m) = -computeFeasibilityMargin(x, y_b, z  - delta, pitch, psi);

//       // std::cout << "y.elem(m) at index " << m << "is " << y.elem(m) << "\n";

//       // std::cout << "margin.states.feetPos[0] " << margin.states.feetPos[0] << "\n";
//       // std::cout << "margin.states.feetPos[1] " << margin.states.feetPos[1] << "\n";
//       // std::cout << "margin.states.feetPos[2] " << margin.states.feetPos[2] << "\n";
//       // std::cout << "margin.states.feetPos[3] " << margin.states.feetPos[3] << "\n";

//       // std::cout << " margin.states.gravityRot " << margin.states.gravityRot << "\n";
//       // std::cout << "pitch " << pitch << "\n";
//       // std::cout << "psi " << psi << "\n";
//       // std::cout << "x " << x << "\n";
//       // std::cout << "y " << y << "\n";
//       // std::cout << "z " << z << "\n";

//       Transform inputTransform;

//       /// Prepare input in right frames
//       // Compute rotation matrix (from base to world) from roll/pitch/yaw - Intrinsic rotations.
//       Eigen::Matrix3d w_R_b = inputTransform.eulerToRotationMatrix(0.0, pitch, psi); // Make function of base rotation, pitch and yaw
//       Eigen::Matrix3d b_R_w = w_R_b.transpose();
//       margin.states.gravityRot = w_R_b.row(2);

//       // Convert to the base frame
//       margin.states.feetPos.resize(4);
//       // margin.states.feetPos[0] = b_R_w * (Eigen::Vector3d(0.44, 0.34, 0.0) - Eigen::Vector3d(x, y_b, z));   // LF Foot Position
//       // margin.states.feetPos[1] = b_R_w * (Eigen::Vector3d(0.44, -0.34, 0.0) - Eigen::Vector3d(x, y_b, z));  // RF Foot Position
//       // margin.states.feetPos[2] = b_R_w * (Eigen::Vector3d(-0.44, 0.34, 0.0) - Eigen::Vector3d(x, y_b, z));  // LH Foot Position
//       // margin.states.feetPos[3] = b_R_w * (Eigen::Vector3d(-0.44, -0.34, 0.0) - Eigen::Vector3d(x, y_b, z)); // RH Foot Position.

//       // margin.states.feetPos[0] = b_R_w * (Eigen::Vector3d(0.55, 0.34, 0.0) - Eigen::Vector3d(x, y_b, z));   // LF Foot Position
//       // margin.states.feetPos[1] = b_R_w * (Eigen::Vector3d(0.54, -0.32, 0.0) - Eigen::Vector3d(x, y_b, z));  // RF Foot Position
//       // margin.states.feetPos[2] = b_R_w * (Eigen::Vector3d(-0.33, 0.30, 0.0) - Eigen::Vector3d(x, y_b, z));  // LH Foot Position
//       // margin.states.feetPos[3] = b_R_w * (Eigen::Vector3d(-0.41, -0.4, 0.0) - Eigen::Vector3d(x, y_b, z)); // RH Foot Position.

//       margin.states.feetPos[0] = b_R_w * (Eigen::Vector3d(0.367, 0.32, 0.0) - Eigen::Vector3d(x, y_b, z));   // LF Foot Position
//       margin.states.feetPos[1] = b_R_w * (Eigen::Vector3d(0.47, -0.41, 0.0) - Eigen::Vector3d(x, y_b, z));  // RF Foot Position
//       margin.states.feetPos[2] = b_R_w * (Eigen::Vector3d(-0.44, 0.37, 0.0) - Eigen::Vector3d(x, y_b, z));  // LH Foot Position
//       margin.states.feetPos[3] = b_R_w * (Eigen::Vector3d(-0.41, -0.36, 0.0) - Eigen::Vector3d(x, y_b, z)); // RH Foot Position.

//       // margin.states.feetPos[0] = b_R_w * (Eigen::Vector3d(0.38, 0.34, 0.2) - Eigen::Vector3d(x, y_b, z));   // LF Foot Position
//       // margin.states.feetPos[1] = b_R_w * (Eigen::Vector3d(0.46, -0.42, 0.2) - Eigen::Vector3d(x, y_b, z));  // RF Foot Position
//       // margin.states.feetPos[2] = b_R_w * (Eigen::Vector3d(-0.44, 0.37, 0.0) - Eigen::Vector3d(x, y_b, z));  // LH Foot Position
//       // margin.states.feetPos[3] = b_R_w * (Eigen::Vector3d(-0.41, -0.36, 0.0) - Eigen::Vector3d(x, y_b, z)); // RH Foot Position.
//       // std::cout << "margin.getFeasibilityGradients(margin.states)" << margin.getFeasibilityGradients(margin.states) << "\n";
//       Eigen::MatrixXd footGradient = margin.getFeasibilityGradients(margin.states);

//       // if (!!J)
//       // {
//       // Partial derivative wrt to X
//       // Transform inputTransform;

//       /// Prepare input in right frames
//       // Compute rotation matrix (from base to world) from roll/pitch/yaw - Intrinsic rotations.
//       // Eigen::Matrix3d w_R_b = inputTransform.eulerToRotationMatrix(0.0, pitch, psi); // Make function of base rotation, pitch and yaw
//       // Eigen::Matrix3d b_R_w = w_R_b.transpose();
//       //  J.elem(m, i - 1) = - (marginPlusZ - marginMinusZ)/(2*delta) ;
//       Eigen::Vector3d dMdG(footGradient(0), footGradient(1), footGradient(2));
//       Eigen::Vector3d dGdTheta(-cos(psi) * cos(pitch), cos(pitch) * sin(psi), -sin(pitch));
//       Eigen::Vector3d dGdPsi(sin(psi) * sin(pitch), sin(pitch) * cos(psi), 0.0);

//       // Eigen::Vector3d dMdG(footGradient(0), footGradient(1), footGradient(2));
//       // Eigen::Vector3d dGdTheta(-cos(pitch), cos(pitch), -sin(pitch));
//       // Eigen::Vector3d dGdPsi(cos(psi) + sin(psi), -sin(psi) + cos(psi), 0.0);

//       J.elem(m, i - 4) = (footGradient(3) + footGradient(6) + footGradient(9) + footGradient(12)) * b_R_w(0, 0);
//       std::cout << "GradX " <<  (footGradient(3) + footGradient(6) + footGradient(9) + footGradient(12)) * b_R_w(0, 0) << "\n";
//       std::cout << "X " <<  x << "\n";

//       m++;

//       J.elem(m, i - 3) = (footGradient(4) + footGradient(7) + footGradient(10) + footGradient(13)) * b_R_w(1, 1);
//       m++;

//       J.elem(m, i - 2) = -dMdG.dot(dGdPsi);
//       m++;

//       J.elem(m, i - 1) = (footGradient(5) + footGradient(8) + footGradient(11) + footGradient(14)) * b_R_w(2, 2);
//       m++;

//       // J.elem(m, i - 4) = (footGradient(15) + footGradient(18) + footGradient(21) + footGradient(24)) * b_R_w(0, 0);
//       // m++;
//       // J.elem(m, i - 3) = (footGradient(16) + footGradient(19) + footGradient(22) + footGradient(25)) * b_R_w(1, 1);
//       // m++;
//       // J.elem(m, i - 2) = dMdG.dot(dGdPsi);
//       // m++;
//       // J.elem(m, i - 1) = (footGradient(17) + footGradient(20) + footGradient(23) + footGradient(26)) * b_R_w(2, 2);
//       // m++;
//       // J.elem(m, i) = 1000000 * (marginPlusPitch - marginMinusPitch);
//       // m++;

//       J.elem(m, i) = -dMdG.dot(dGdTheta);
//       m++;

//       // J.elem(m, i - 1) = 0.0;
//       // }
//       // J.elem(m, i) = 2.0 * qi + 0.4;
//       // m++;
//     }

//     // std::cout << "J " << J << "\n";
//   }
// }

// double F_feasibilityMarginConstraint::computeFeasibilityMarginConstraint(double x, double y, double z, double pitch, double psi)
// {
//   Transform inputTransform;

//   /// Prepare input in right frames
//   // Compute rotation matrix (from base to world) from roll/pitch/yaw - Intrinsic rotations.
//   Eigen::Matrix3d w_R_b = inputTransform.eulerToRotationMatrix(0.0, pitch, psi); // Make function of base rotation, pitch and yaw
//   Eigen::Matrix3d b_R_w = w_R_b.transpose();
//   margin.states.gravityRot = w_R_b.row(2);

//   // Convert to the base frame
//   margin.states.feetPos.resize(4);
//   // margin.states.feetPos[0] = b_R_w * (Eigen::Vector3d(0.44, 0.34, 0.0) - Eigen::Vector3d(x, y, z));   // LF Foot Position
//   // margin.states.feetPos[1] = b_R_w * (Eigen::Vector3d(0.44, -0.34, 0.0) - Eigen::Vector3d(x, y, z));  // RF Foot Position
//   // margin.states.feetPos[2] = b_R_w * (Eigen::Vector3d(-0.44, 0.34, 0.0) - Eigen::Vector3d(x, y, z));  // LH Foot Position
//   // margin.states.feetPos[3] = b_R_w * (Eigen::Vector3d(-0.44, -0.34, 0.0) - Eigen::Vector3d(x, y, z)); // RH Foot Position.

//       // margin.states.feetPos[0] = b_R_w * (Eigen::Vector3d(0.55, 0.34, 0.0) - Eigen::Vector3d(x, y, z));   // LF Foot Position
//       // margin.states.feetPos[1] = b_R_w * (Eigen::Vector3d(0.54, -0.32, 0.0) - Eigen::Vector3d(x, y, z));  // RF Foot Position
//       // margin.states.feetPos[2] = b_R_w * (Eigen::Vector3d(-0.33, 0.30, 0.0) - Eigen::Vector3d(x, y, z));  // LH Foot Position
//       // margin.states.feetPos[3] = b_R_w * (Eigen::Vector3d(-0.41, -0.4, 0.0) - Eigen::Vector3d(x, y, z)); // RH Foot Position.



//       margin.states.feetPos[0] = b_R_w * (Eigen::Vector3d(0.367, 0.32, 0.0) - Eigen::Vector3d(x, y, z));   // LF Foot Position
//       margin.states.feetPos[1] = b_R_w * (Eigen::Vector3d(0.47, -0.41, 0.0) - Eigen::Vector3d(x, y, z));  // RF Foot Position
//       margin.states.feetPos[2] = b_R_w * (Eigen::Vector3d(-0.44, 0.37, 0.0) - Eigen::Vector3d(x, y, z));  // LH Foot Position
//       margin.states.feetPos[3] = b_R_w * (Eigen::Vector3d(-0.41, -0.36, 0.0) - Eigen::Vector3d(x, y, z)); // RH Foot Position.

//   return margin.getFeasibilityMargin(margin.states);
// }

// void F_feasibilityMarginConstraint::init()
// {
//   std::string NewtorkPath = "/home/simple/ros_ws_lgp/src/reachability-aware-LGP/tamp_lgp/../../networks_minimal";
//   std::string parametersDirectory(NewtorkPath + "/examples/parameters/");
//   std::cout << "INIT_FEASIBILITY" << "\n";

//   margin.loadDataScalersFromFile(parametersDirectory + "hyqreal/dataScalers.txt");
// }

//===========================================================================

uintA getNonSwitchedFrames(const FrameL &A, const FrameL &B)
{
  uintA nonSwitchedFrames;
  CHECK_EQ(A.N, B.N, "");

  for (uint i = 0; i < A.N; i++)
  {
    rai::Frame *f0 = A.elem(i);
    rai::Frame *f1 = B.elem(i);
    if (!f0->joint || !f1->joint)
      continue;
    if (f0->joint->type != f1->joint->type)
      continue;
    if (f0->joint->mimic || f1->joint->mimic)
      continue;
    if (f0->ID - f0->parent->ID != f1->ID - f1->parent->ID)
      continue; // comparing the DIFFERENCE in IDs between parent and joint
    nonSwitchedFrames.append(i);
  }
  return nonSwitchedFrames;
}
