/*  ------------------------------------------------------------------
    Copyright (c) 2011-2020 Marc Toussaint
    email: toussaint@tu-berlin.de

    This code is distributed under the MIT License.
    Please see <root-path>/LICENSE for details.
    --------------------------------------------------------------  */

#include "../Logic/relationalMachine.h"
#include "../Core/thread.h"

struct ServiceRAP {
  unique_ptr<struct sServiceRAP> self;

  ServiceRAP();
  ~ServiceRAP();
};
