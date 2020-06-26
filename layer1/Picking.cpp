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
 * Value to put in the "unused bits" part of a color channel.
 * May serve several purposes:
 * - validation (check bits)
 * - no-pick vs. through-pick
 * - rounding
 *
 * Disclaimer: I'm not sure if rounding is ever needed and if rounding and
 * validation should be mutually exclusive.
 */
constexpr unsigned make_check_value(unsigned bits)
{
  return 0x80u >> bits;
}

/**
 * Return true if this RGBA color looks valid, and false if it should
 * be ignored.
 */
bool PickColorConverter::validateCheckBits(const channel_t* rgba) const
{
  for (unsigned i = 0; i != 4; ++i) {
    assert(m_rgba_and_check_bits[i] >= m_rgba_bits[i]);

    // mask to stamp out the validation bits
    const channel_t check_mask = (0xFFu >> m_rgba_bits[i]) & //
                                 ~(0xFFu >> m_rgba_and_check_bits[i]);

    const channel_t check_value = make_check_value(m_rgba_bits[i]);

    if ((rgba[i] & check_mask) != (check_value & check_mask)) {
      // found antialiased bit (or cPickableNoPick/cPickableThrough which can
      // also be ignored)
      return false;
    }
  }

  return true;
}

/**
 * Set the number of bits per picking channel
 * @param rgba_bits Number of available bits per channel
 * @param max_check_bits Maximum number of bits (per channel) to use for
 * validation
 */
void PickColorConverter::setRgbaBits(const int* rgba_bits, int max_check_bits)
{
  for (unsigned i = 0; i != 4; ++i) {
    m_rgba_bits[i] = std::min<unsigned>(MAX_BITS, rgba_bits[i]);

    // Use at most half of the pixels for validation
    int const check_bits = std::min(m_rgba_bits[i] / 2, max_check_bits);

    m_rgba_and_check_bits[i] = m_rgba_bits[i];
    m_rgba_bits[i] = std::max(0, int(m_rgba_bits[i]) - check_bits);
  }

  // Use one bit to distinguish no-pick and through-pick.
  // This bit only has to reach the shader, it's not relevant for the
  // output buffer. This means it can be part of the check bits, but
  // doesn't have to be. It can also be beyond the color depth of the
  // output buffer (In theory, I have no suitable system to test this
  // assumption).
  m_rgba_bits[3] = std::min<unsigned>(MAX_BITS - 1, m_rgba_bits[3]);
}

/**
 * Get the picking index from color
 * @param rgba RGBA color
 * @return Picking index
 */
PickColorConverter::index_t PickColorConverter::indexFromColor(
    const channel_t* rgba) const
{
  if (!validateCheckBits(rgba)) {
    return 0;
  }

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

    rgba[i] |= make_check_value(m_rgba_bits[i]);

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
  rgba[3] = make_check_value(m_rgba_bits[3]);
  assert(rgba[3] != 0);
}

/**
 * Get pick-through color
 * @param[out] rgba RGBA color
 *
 * Note: This only works if the shader discards full-transparent pixels (e.g. may not
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
