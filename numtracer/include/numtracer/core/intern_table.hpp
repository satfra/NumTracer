/// @file core/intern_table.hpp
/// @brief Open-addressed interning index: maps a key to its position in a caller-owned vector.
///
/// Used by @ref numtracer::network::GlobalEnv (`codegen/gen.hpp`, the shared `f[]` symbol table) and
/// @ref numtracer::network::rdetail::RBuilder (`codegen/real_cse.hpp`, value numbering of the real
/// SSA). Both own the key vector themselves, because the rest of the engine reads it directly
/// (`syms`, `ins`); this index only adds the hash lookup on top. Keys must therefore be appended
/// only through @ref InternTable::intern, or the index goes stale.
///
/// Ids are positions in that vector, so they are assigned in first-seen order and do not depend on
/// the hash at all.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace numtracer
{

  /// @brief Linear-probing hash index over a `std::vector<Key>`, kept at load factor < 0.7.
  ///
  /// `Hash` maps a key to `std::uint64_t`, `Eq` compares two keys; both are stateless.
  template <class Key, class Hash, class Eq> struct InternTable {
    std::vector<int> bucket; ///< Open-addressed index: position in the key vector, or -1 if empty.
    std::size_t mask = 0;    ///< `bucket.size()-1` (power of two); 0 while empty.

    /// @brief The position of @p key in @p keys, appending it on a miss.
    constexpr int intern(std::vector<Key> &keys, const Key &key)
    {
      const std::uint64_t h = Hash{}(key);
      if (mask) {
        std::size_t p = h & mask;
        while (bucket[p] != -1) {
          if (Eq{}(keys[bucket[p]], key)) return bucket[p];
          p = (p + 1) & mask;
        }
      }
      if ((keys.size() + 1) * 10 >= (mask + 1) * 7) rehash(keys, mask == 0 ? 16 : (mask + 1) * 2);
      const int s = static_cast<int>(keys.size());
      keys.push_back(key);
      insert_slot(h, s);
      return s;
    }

  private:
    constexpr void insert_slot(std::uint64_t h, int s)
    {
      std::size_t p = h & mask;
      while (bucket[p] != -1)
        p = (p + 1) & mask;
      bucket[p] = s;
    }
    constexpr void rehash(const std::vector<Key> &keys, std::size_t cap)
    {
      bucket.assign(cap, -1);
      mask = cap - 1;
      for (std::size_t s = 0; s < keys.size(); ++s)
        insert_slot(Hash{}(keys[s]), static_cast<int>(s));
    }
  };

} // namespace numtracer
