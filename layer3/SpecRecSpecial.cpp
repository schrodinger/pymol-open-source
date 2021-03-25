#include "SpecRecSpecial.h"
#include "SpecRec.h"

OrderRec::OrderRec(std::string recName_, std::size_t pos_)
  : name(std::move(recName_))
  , pos(pos_)
{}

DiscardedRec::DiscardedRec(SpecRec* rec_, std::size_t pos_)
  : rec(rec_)
  , pos(pos_)
{}
