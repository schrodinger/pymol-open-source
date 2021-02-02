/*
 * Color picking
 *
 * (C) Schrodinger, Inc.
 */

#pragma once

#include <vector>

//! With shaders we have to pre-allocate picking color buffers, so
//! limit this to 2 passes.
constexpr auto SHADER_PICKING_PASSES_MAX = 2;

//! Transparency cutoff for cPickableThrough, smaller values are considered
//! opaque.
constexpr auto PICKABLE_THROUGH_CUTOFF = 0.1f;

//! Valid values for the `transparency_picking_mode` setting
enum {
  cTransparencyPickingModePickable = 1,
  cTransparencyPickingModeAuto = 2, // 0 has same effect
};

enum cPickable_t {
  cPickableAtom = -1,
  cPickableLabel = -2,
  cPickableGadget = -3,
  cPickableNoPick = -4,
  cPickableThrough = -5,
};

/**
 * Identifier for a pickable item like an atom, half-bond or label, within a
 * pick context (object state).
 */
struct Pickable {
  //! Primary index (e.g. atom index within object)
  unsigned int index;

  //! Secondary index (>= 0, e.g. bond index) or special value (::cPickable_t)
  int bond;

  // comparison
  bool operator==(const Pickable& rhs) const
  {
    return index == rhs.index && bond == rhs.bond;
  }
};

namespace pymol
{
struct CObject;
}

/**
 * Generic object state identifier
 */
struct PickContext {
  pymol::CObject* object = nullptr;
  int state;

  // comparison
  bool operator==(const PickContext& rhs) const
  {
    return object == rhs.object && state == rhs.state;
  }
};

/**
 * Global identifier for a pickable item like an atom, half-bond or label
 */
struct Picking {
  Pickable src;
  PickContext context;

  // comparison
  bool operator==(const Picking& rhs) const
  {
    return src == rhs.src && context == rhs.context;
  }
  bool operator!=(const Picking& rhs) const { return !(*this == rhs); }
};

/**
 * For picking everything within a rectangle (Shift+click drag selection)
 */
struct Multipick {
  int x, y, w, h;
  std::vector<Picking> picked;
};

/**
 * Color picking index conversions
 */
class PickColorConverter
{
  typedef unsigned char channel_t;
  typedef unsigned int index_t;

  // 12 bit default
  unsigned char m_rgba_bits[4]{4, 4, 4, 0};

  // Number of used bits plus check bits
  unsigned char m_rgba_and_check_bits[4]{4, 4, 4, 0};

public:
  // (only public for testing)
  bool validateCheckBits(const channel_t* rgba) const;

  //! Set the number of bits for each channel
  void setRgbaBits(const int* rgba_bits, int max_check_bits = 0);

  //! Get the sum of bits from all channels
  unsigned getTotalBits() const
  {
    return m_rgba_bits[0] + m_rgba_bits[1] + m_rgba_bits[2] + m_rgba_bits[3];
  }

  //! Picking index from picked color
  index_t indexFromColor(const channel_t* rgba) const;

  //! RGBA color from picking index
  void colorFromIndex(channel_t* rgba, index_t idx) const;

  //! RGBA color for "no-pick" (blocks pick, but doesn't pick anything)
  void colorNoPick(channel_t* rgba) const;

  //! RGBA color for "pick-through" (doesn't block, can pick items below)
  void colorPickThrough(channel_t* rgba) const;
};

class PickColorManager : public PickColorConverter
{
  unsigned m_count = 0;
  std::vector<Picking> m_identifiers;

public:
  //! picking pass
  unsigned m_pass = 0;

  //! whether `m_identifiers` and `m_count` are complete
  bool m_valid = false;

  //! Running index while `!m_valid`, equal to `m_identifiers.size()` at the end
  //! of a picking pass or if `m_valid`.
  unsigned count() const { return m_count; }

  //! Reset the count, but don't invalidate identifiers (they will be identical
  //! in the next pass)
  void resetCount()
  {
    m_count = 0;
    m_valid = false;
  }

  //! Invalidate identifiers
  void invalidate()
  {
    if (!m_valid) {
      // protect against calls during rendering
      return;
    }
    m_identifiers.clear();
    m_valid = false;
  }

  bool pickColorsValid() const { return m_valid; }

  void colorNext(unsigned char* color, const PickContext* context,
      unsigned int index, int bond);

  const Picking* getIdentifier(unsigned index) const;
};
