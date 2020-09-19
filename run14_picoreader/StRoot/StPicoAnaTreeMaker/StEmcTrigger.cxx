#include "StEmcTrigger.h"
#include "StPicoDstMaker/StPicoConstants.h"
#include "StMessMgr.h"
#include <climits>

ClassImp(StEmcTrigger)

//----------------------------------------------------------------------------------
StEmcTrigger::StEmcTrigger()
{
  Clear();
}

//----------------------------------------------------------------------------------
StEmcTrigger::StEmcTrigger(int flag, int id, int adc, int eId, int adc0)
{
  Clear();

  if(flag<0) mFlag = 0;
  if(id  <0) mId   = 0;
  if(adc <0) mAdc  = 0;
  if(eId <0) mEId  = USHRT_MAX;
  if(adc0<0) mAdc0 = 0;

  mFlag = (flag>UCHAR_MAX)  ? UCHAR_MAX  : (UChar_t)flag;
  mId   = (id  >USHRT_MAX) ? USHRT_MAX : (UShort_t)id;
  mAdc  = (adc >USHRT_MAX) ? USHRT_MAX : (UShort_t)adc;
  mEId = (eId >USHRT_MAX) ? USHRT_MAX : (UShort_t)eId;
  mAdc0  = (adc0 >USHRT_MAX) ? USHRT_MAX : (UShort_t)adc0;
}

//----------------------------------------------------------------------------------
StEmcTrigger::~StEmcTrigger()
{ /* noop */ }

//----------------------------------------------------------------------------------
void StEmcTrigger::Clear(const Option_t* opt)
{
  mFlag = 0;
  mId = 0;
  mAdc = 0;
  mEId = 0;
  mAdc0 = 0;

}
//----------------------------------------------------------------------------------
void StEmcTrigger::Print(const Char_t *option) const {
  LOG_INFO << " Flag = " << mFlag << " Id = " << mId << " Adc = " << mAdc << " mEId = "<<mEId<< endm;
}
