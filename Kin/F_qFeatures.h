/*  ------------------------------------------------------------------
    Copyright (c) 2011-2020 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#pragma once

#include "feature.h"

//===========================================================================

struct F_qItself : Feature
{
  enum PickMode
  {
    byJointNames,
    byFrameNames,
    byExcludeJointNames,
    allActiveJoints
  };

  bool relative_q0; ///< if true, absolute values are given relative to Joint::q0

  F_qItself(bool relative_q0 = false);
  F_qItself(PickMode pickMode, const StringA &picks, const rai::Configuration &C, bool relative_q0 = false);
  F_qItself(const uintA &_selectedFrames, bool relative_q0 = false);

  virtual void phi(arr &y, arr &J, const rai::Configuration &C);
  virtual void phi2(arr &y, arr &J, const FrameL &F);
  virtual uint dim_phi(const rai::Configuration &C);
  virtual uint dim_phi2(const FrameL &F);

private:
  std::map<rai::Configuration *, uint> dimPhi;
};

//===========================================================================

struct F_qZeroVel : Feature
{
  virtual void phi2(arr &y, arr &J, const FrameL &F);
  virtual uint dim_phi2(const FrameL &F);
};

//===========================================================================

struct F_qLimits2 : Feature
{
  virtual void phi2(arr &y, arr &J, const FrameL &F);
  virtual uint dim_phi2(const FrameL &F);
};

//===========================================================================

struct F_qLimits : Feature
{
  virtual void phi2(arr &y, arr &J, const FrameL &F);
  virtual uint dim_phi2(const FrameL &F);
};

//===========================================================================

struct F_qQuaternionNorms : Feature
{
  virtual void phi2(arr &y, arr &J, const FrameL &F);
  virtual uint dim_phi2(const FrameL &F);

  void setAllActiveQuats(const rai::Configuration &C);
};

//===========================================================================

struct F_qTime : Feature
{
  virtual void phi2(arr &y, arr &J, const FrameL &F);
  virtual uint dim_phi2(const FrameL &F) { return 1; }
};

//===========================================================================

// struct F_parabola : Feature
// {

//   // New method to compute the parabola
//   // virtual void phi2(arr& y, arr& J, const rai::Configuration& C);
//   virtual void phi2(arr &y, arr &J, const FrameL &F);
//   virtual uint dim_phi2(const FrameL &F) { return 5; }
// };

//===========================================================================

// struct F_feasibilityMargin : Feature
// {
//   FeasibilityMargin margin{"/home/simple/ros_ws_lgp/src/networks_minimal/examples/parameters/hyqreal",
//                              {15, 256, 256, 128, 1},
//                              activation.relu};
                             
//   F_feasibilityMargin()
//   {
   
//     init();
//   };
//   void init();
//   double computeFeasibilityMargin(double x, double y, double z, double pitch, double psi);
//   virtual void phi2(arr &y, arr &J, const FrameL &F);
//   virtual uint dim_phi2(const FrameL &F) { return 1; }
// };

// //===========================================================================

// struct F_feasibilityMarginConstraint : Feature
// {
//   FeasibilityMargin margin{"/home/simple/ros_ws_lgp/src/networks_minimal/examples/parameters/hyqreal",
//                              {15, 256, 256, 128, 1},
//                              activation.relu};
                             
//   F_feasibilityMarginConstraint()
//   {
   
//     init();
//   };
//   void init();
//   double computeFeasibilityMarginConstraint(double x, double y, double z, double pitch, double psi);
//   virtual void phi2(arr &y, arr &J, const FrameL &F);
//   virtual uint dim_phi2(const FrameL &F) { return 5; }
// };

//===========================================================================

rai::Array<rai::Joint *> getMatchingJoints(const ConfigurationL &Ktuple, bool zeroVelJointsOnly);
uintA getNonSwitchedFrames(const FrameL &A, const FrameL &B);