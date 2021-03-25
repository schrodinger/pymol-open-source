#pragma once

#include <string>

class SpecRec;

struct OrderRec {
  std::string name;
  std::size_t pos;
  OrderRec(std::string rec_, std::size_t pos_);
};

/**
 * Contains a discarded rec and its previous position
 * on the panel list
 *
 * @note if retrived from a function (ExecutiveDelete w/ save),
 * the caller is responsible for calling ExecutivePurgeSpecRec
 * w/o saving on the contained rec.
 */

struct DiscardedRec {
  SpecRec* rec;
  std::size_t pos;
  DiscardedRec(SpecRec* rec, std::size_t pos);
};
