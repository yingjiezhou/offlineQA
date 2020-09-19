#include <iostream>
#include <bitset>

#include "StTriggerData.h"
#include "StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuDst.h" 
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuMtdHeader.h"
#include "StPicoDstMaker/StPicoMtdTrigger.h"
#include "StMtdTrigger.h"

ClassImp(StMtdTrigger)

//----------------------------------------------------------------------------------
StMtdTrigger::StMtdTrigger()
{
  memset(mQTtacSum,0,sizeof(mQTtacSum));
  memset(mMT101Tac,0,sizeof(mMT101Tac));
  memset(mMT101Id,0,sizeof(mMT101Id));
  mTF201TriggerBit = 0;
}

//----------------------------------------------------------------------------------
StMtdTrigger::StMtdTrigger(StPicoMtdTrigger *mtd)
{

	for(int i=0;i<4;i++){
		for(int j=0;j<8;j++){
			mQTtacSum[i][j] = mtd->getQTtacSum(i+1,j+1);
		}
		for(int j=0;j<2;j++){
			mMT101Tac[i][j] = mtd->getMT101Tac(i+1,j);
			mMT101Id[i][j] = mtd->getMT101Id(i+1,j);
		}
	}

	mTF201TriggerBit = mtd->getTF201TriggerBit();
//    mShouldHaveRejectEvent = mtd->shouldHaveRejectEvent();
}

//----------------------------------------------------------------------------------
StMtdTrigger::~StMtdTrigger()
{
}


//----------------------------------------------------------------------------------
void StMtdTrigger::getMaximumQTtac(const Int_t qt, Int_t& pos1, Int_t& pos2)
{
  pos1 = 0;
  pos2 = 0;

  UShort_t max1 = 0, max2 = 0;

  for(Int_t i=0; i<8; i++)
    {
      if(max1 < mQTtacSum[qt-1][i])
	{
	  max2 = max1;
	  pos2 = pos1;
	  max1 = mQTtacSum[qt-1][i];
	  pos1 = i+1;
	}
      else if(max2 < mQTtacSum[qt-1][i])
	{
	  max2 = mQTtacSum[qt-1][i];
	  pos2 = i+1;
	}
    }
}
