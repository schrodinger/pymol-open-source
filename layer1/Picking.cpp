/*
 * Color picking
 *
 * (C) Schrodinger, Inc.
 */

#include "Picking.h"

#include <algorithm>
#include <cassert>

constexpr auto MAX_BITS = 8u;

/**
 * Set the number of bits per picking channel
 */
void PickColorConverter::setRgbaBits(const int* rgba_bits)
{
  m_rgba_bits[0] = std::min<unsigned>(MAX_BITS, rgba_bits[0]);
  m_rgba_bits[1] = std::min<unsigned>(MAX_BITS, rgba_bits[1]);
  m_rgba_bits[2] = std::min<unsigned>(MAX_BITS, rgba_bits[2]);

  // we use one bit to distinguish no-pick and through-pick
  m_rgba_bits[3] = std::min<unsigned>(MAX_BITS - 1, rgba_bits[3]);
}

/**
 * Get the picking index from color
 * @param rgba RGBA color
 * @return Picking index
 */
PickColorConverter::index_t PickColorConverter::indexFromColor(
    const channel_t* rgba) const
{
  unsigned idx = 0;
  unsigned bits = 0;
  for (unsigned i = 0; i != 4; ++i) {
    idx |= (unsigned(rgba[i]) >> (MAX_BITS - m_rgba_bits[i])) << bits;
    bits += m_rgba_bits[i];
  }
  return idx;
}

/**
 * Get color from picking index
 * @param[out] rgba RGBA color
 * @param idx Picking index
 */
void PickColorConverter::colorFromIndex(channel_t* rgba, index_t idx) const
{
  unsigned bits = 0;
  for (unsigned i = 0; i != 4; ++i) {
    rgba[i] = ((idx >> bits) & 0xFFu) << (MAX_BITS - m_rgba_bits[i]);

    // "rounding" by setting the next bit (like adding 0.5 - values are
    // effectively floored in indexFromColor()) - may not be necessary?
    rgba[i] |= 0x80u >> m_rgba_bits[i];

    bits += m_rgba_bits[i];
  }
}

/**
 * Get no-pick color
 * @param[out] rgba RGBA color
 */
void PickColorConverter::colorNoPick(channel_t* rgba) const
{
  rgba[0] = 0;
  rgba[1] = 0;
  rgba[2] = 0;
  rgba[3] = 0x80u >> m_rgba_bits[3]; // lowest bit
}

/**
 * Get pick-through color
 * @param[out] rgba RGBA color
 *
 * Note: This only works if the shader discards full-opaque pixels (e.g. may not
 * work in immediate mode)
 */
void PickColorConverter::colorPickThrough(channel_t* rgba) const
{
  rgba[0] = 0;
  rgba[1] = 0;
  rgba[2] = 0;
  rgba[3] = 0;
}

/**
 * Get the color for (context, index, bond, current picking pass).
 */
void PickColorManager::colorNext(unsigned char* color,
    const PickContext* context, unsigned int index, int bond)
{
  if (bond == cPickableNoPick) {
    colorNoPick(color);
    return;
  }

  if (bond == cPickableThrough) {
    // TODO has no effect on immediate mode (fragment not discarded)
    colorPickThrough(color);
    return;
  }

  const Picking p_new = {{index, bond}, *context};

  assert(m_count <= m_identifiers.size());

  // compare with previous pick color, increment if different
  if (m_count == 0 || m_identifiers[m_count - 1] != p_new) {
    ++m_count;
  }

  unsigned j = m_count;

  if (m_pass > 0) {
    assert(m_count <= m_identifiers.size());
    j >>= getTotalBits() * m_pass;
  } else if (m_count == m_identifiers.size() + 1) {
    m_identifiers.push_back(p_new);
  }

  // if this assertion fails, then invalidate() was not called, but should be
  assert(m_identifiers[m_count - 1] == p_new);

  colorFromIndex(color, j);
}

/**
 * Get identifier for the 1-based picking color index.
 * Returns NULL if the index is out of bounds.
 */
const Picking* PickColorManager::getIdentifier(unsigned index) const
{
  if (index > 0 && index <= m_identifiers.size()) {
    return m_identifiers.data() + index - 1;
  }
  return nullptr;
}
