/// @file expr.hpp
/// @brief Compile-time expression algebra with structural zero/one pruning.
///
/// An *expression type* is an empty C++ type carrying no runtime state; the
/// type itself *is* the expression. Building an expression composes these
/// types, and the normalizing combinators @ref numtracer::expr::add / @ref
/// numtracer::expr::mul apply the algebra's identity, annihilator and
/// constant-folding laws **at compile time**. The effect is *structural
/// pruning*: a branch that touches a structural @ref numtracer::expr::Zero (for
/// example a zero entry of a gamma matrix) disappears from the resulting type
/// entirely. After a contraction, the surviving type for each output slot is a
/// pruned sum-of-products of runtime leaves (@ref numtracer::expr::Var), which
/// @ref numtracer::expr::evaluate turns into a number given the per-grid-point
/// values.
///
/// Node types:
///   - @ref numtracer::expr::Zero — additive identity, multiplicative
///   annihilator
///     (evaluates to 0).
///   - @ref numtracer::expr::One — multiplicative identity (evaluates to 1).
///   - @ref numtracer::expr::Lit — any other compile-time complex constant.
///   - @ref numtracer::expr::Var — a runtime value, looked up as `env[Id]`.
///   - @ref numtracer::expr::Add, @ref numtracer::expr::Mul — binary nodes.
///
/// Rewrite laws applied by `add`/`mul`: `Zero*x = Zero`, `One*x = x`,
/// `Zero+x = x`, and constant folding of `Lit*Lit` / `const+const`. Folding
/// re-normalizes back to `Zero`/`One` when the constant is 0/1, so e.g. `1*x`
/// reduces to `x` and the identities keep firing downstream.
#pragma once

#include "numtracer/core/cx.hpp" // Cx (constexpr complex)

#include <type_traits>

namespace numtracer::expr {

// ---- node types (all empty) ------------------------------------------------

/// @brief The additive identity and multiplicative annihilator (evaluates to
/// 0).
struct Zero {};
/// @brief The multiplicative identity (evaluates to 1).
struct One {};
/// @brief A compile-time complex constant other than 0 or 1.
/// @tparam C The constant value.
template <Cx C> struct Lit {};
/// @brief A runtime leaf, evaluated as `env[Id]`.
/// @tparam Id Index into the runtime environment array (see the "env id
/// registry"
///         conventions documented on the sector builders that emit `Var`s).
template <int Id> struct Var {};
/// @brief Binary addition node `L + R`.
/// @tparam L Left operand expression type.
/// @tparam R Right operand expression type.
template <class L, class R> struct Add {};
/// @brief Binary multiplication node `L * R`.
/// @tparam L Left operand expression type.
/// @tparam R Right operand expression type.
template <class L, class R> struct Mul {};

// ---- traits ----------------------------------------------------------------

/// @brief `true` if `T` is @ref numtracer::expr::Zero.
/// @tparam T The type to test.
template <class T> inline constexpr bool is_zero = std::is_same_v<T, Zero>;
/// @brief `true` if `T` is @ref numtracer::expr::One.
/// @tparam T The type to test.
template <class T> inline constexpr bool is_one = std::is_same_v<T, One>;

/// @brief Whether an expression type has a known compile-time value, and what
/// it is.
///
/// Specialized so that `const_val<T>::is` is `true` for @ref
/// numtracer::expr::Zero,
/// @ref numtracer::expr::One and @ref numtracer::expr::Lit, in which case
/// `const_val<T>::value` holds the @ref numtracer::Cx value. For all other node
/// types `is` is `false`.
/// @tparam T The expression type to inspect.
template <class T> struct const_val {
  static constexpr bool is =
      false; ///< `true` when the value is known at compile time.
};
/// @brief Constant trait for @ref numtracer::expr::Zero (value `0`).
template <> struct const_val<Zero> {
  static constexpr bool is = true; ///< The value is known at compile time.
  static constexpr Cx value{0, 0}; ///< The constant value.
};
/// @brief Constant trait for @ref numtracer::expr::One (value `1`).
template <> struct const_val<One> {
  static constexpr bool is = true; ///< The value is known at compile time.
  static constexpr Cx value{1, 0}; ///< The constant value.
};
/// @brief Constant trait for @ref numtracer::expr::Lit (value `C`).
/// @tparam C The literal value.
template <Cx C> struct const_val<Lit<C>> {
  static constexpr bool is = true; ///< The value is known at compile time.
  static constexpr Cx value = C;   ///< The constant value.
};

/// @brief Marks the expression node types so the operators below apply only to
/// them.
///
/// `false` by default; specialized to `true` for each of the six node types.
/// This is what stops `operator*` / `operator+` from hijacking ordinary
/// arithmetic (e.g. `Cx * Cx`).
/// @tparam T The type to test.
template <class T> inline constexpr bool is_expr = false;
template <>
inline constexpr bool is_expr<Zero> =
    true; ///< @ref Zero is an expression node.
template <>
inline constexpr bool is_expr<One> = true; ///< @ref One is an expression node.
template <Cx C>
inline constexpr bool is_expr<Lit<C>> =
    true; ///< @ref Lit is an expression node.
template <int Id>
inline constexpr bool is_expr<Var<Id>> =
    true; ///< @ref Var is an expression node.
template <class L, class R>
inline constexpr bool is_expr<Add<L, R>> =
    true; ///< @ref Add is an expression node.
template <class L, class R>
inline constexpr bool is_expr<Mul<L, R>> =
    true; ///< @ref Mul is an expression node.

/// @brief Number of nodes (leaves plus `Add`/`Mul`) in a pruned expression
/// type.
///
/// A direct measure of how much the structural zero/one pruning shrank the
/// surviving tree. Leaves (`Zero`/`One`/`Lit`/`Var`) count as 1.
/// @tparam E The expression type to measure.
template <class E>
inline constexpr int node_count = 1; // Zero/One/Lit/Var are leaves
template <class L, class R>
inline constexpr int node_count<Add<L, R>> =
    1 + node_count<L> + node_count<R>; ///< @ref Add node, 1 plus both children.
template <class L, class R>
inline constexpr int node_count<Mul<L, R>> =
    1 + node_count<L> + node_count<R>; ///< @ref Mul node, 1 plus both children.

/// @brief Normalize a compile-time complex constant to its canonical node type.
/// @tparam C The constant value.
/// @return `Zero{}` if `C == 0`, `One{}` if `C == 1`, otherwise `Lit<C>{}`.
/// This
///         keeps the identity/annihilator rewrites firing after constant
///         folding.
template <Cx C> constexpr auto lit_norm() {
  if constexpr (C.re == 0.0 && C.im == 0.0)
    return Zero{};
  else if constexpr (C.re == 1.0 && C.im == 0.0)
    return One{};
  else
    return Lit<C>{};
}

// ---- normalizing combinators ----------------------------------------------

/// @brief Multiply two expression types, applying the algebra's rewrite laws.
///
/// Rewrites: `Zero*x = Zero`, `x*Zero = Zero`, `One*x = x`, `x*One = x`, and
/// `Lit*Lit` folds (re-normalized via @ref numtracer::expr::lit_norm).
/// Otherwise yields `Mul<A,B>`.
/// @tparam A Left operand expression type (deduced).
/// @tparam B Right operand expression type (deduced).
/// @return The normalized product expression type.
template <class A, class B> constexpr auto mul(A, B) {
  if constexpr (is_zero<A> || is_zero<B>)
    return Zero{};
  else if constexpr (is_one<A>)
    return B{};
  else if constexpr (is_one<B>)
    return A{};
  else if constexpr (const_val<A>::is && const_val<B>::is) // Lit * Lit -> fold
    return lit_norm<const_val<A>::value * const_val<B>::value>();
  else
    return Mul<A, B>{};
}

/// @brief Add two expression types, applying the algebra's rewrite laws.
///
/// Rewrites: `Zero+x = x`, `x+Zero = x`, and `const+const` folds (re-normalized
/// via @ref numtracer::expr::lit_norm). Otherwise yields `Add<A,B>`.
/// @tparam A Left operand expression type (deduced).
/// @tparam B Right operand expression type (deduced).
/// @return The normalized sum expression type.
template <class A, class B> constexpr auto add(A, B) {
  if constexpr (is_zero<A>)
    return B{};
  else if constexpr (is_zero<B>)
    return A{};
  else if constexpr (const_val<A>::is &&
                     const_val<B>::is) // const + const -> fold
    return lit_norm<const_val<A>::value + const_val<B>::value>();
  else
    return Add<A, B>{};
}

/// @brief `*` sugar for @ref numtracer::expr::mul, restricted to expression
/// nodes.
///
/// The `is_expr` constraint ensures ordinary arithmetic (e.g. `Cx*Cx`) is
/// unaffected.
/// @tparam A Left operand expression type (deduced).
/// @tparam B Right operand expression type (deduced).
/// @param a Left operand (value ignored; types carry the expression).
/// @param b Right operand (value ignored; types carry the expression).
/// @return The normalized product expression type.
template <class A, class B>
  requires(is_expr<A> && is_expr<B>)
constexpr auto operator*(A a, B b) {
  return mul(a, b);
}
/// @brief `+` sugar for @ref numtracer::expr::add, restricted to expression
/// nodes.
/// @tparam A Left operand expression type (deduced).
/// @tparam B Right operand expression type (deduced).
/// @param a Left operand (value ignored; types carry the expression).
/// @param b Right operand (value ignored; types carry the expression).
/// @return The normalized sum expression type.
template <class A, class B>
  requires(is_expr<A> && is_expr<B>)
constexpr auto operator+(A a, B b) {
  return add(a, b);
}

// ---- evaluation: plug runtime values (env[Id]) into the pruned tree --------

/// @brief Evaluate @ref numtracer::expr::Zero.
/// @return `Cx{0,0}`.
constexpr Cx eval(Zero, const Cx *) { return {0, 0}; }
/// @brief Evaluate @ref numtracer::expr::One.
/// @return `Cx{1,0}`.
constexpr Cx eval(One, const Cx *) { return {1, 0}; }
/// @brief Evaluate a @ref numtracer::expr::Lit.
/// @tparam C The literal value.
/// @return `C`.
template <Cx C> constexpr Cx eval(Lit<C>, const Cx *) { return C; }
/// @brief Evaluate a @ref numtracer::expr::Var by environment lookup.
/// @tparam Id The variable's environment index.
/// @param env The runtime environment array.
/// @return `env[Id]`.
template <int Id> constexpr Cx eval(Var<Id>, const Cx *env) { return env[Id]; }
/// @brief Evaluate an @ref numtracer::expr::Add node.
/// @tparam L Left operand expression type.
/// @tparam R Right operand expression type.
/// @param env The runtime environment array.
/// @return The sum of the evaluated children.
template <class L, class R> constexpr Cx eval(Add<L, R>, const Cx *env) {
  return eval(L{}, env) + eval(R{}, env);
}
/// @brief Evaluate a @ref numtracer::expr::Mul node.
/// @tparam L Left operand expression type.
/// @tparam R Right operand expression type.
/// @param env The runtime environment array.
/// @return The product of the evaluated children.
template <class L, class R> constexpr Cx eval(Mul<L, R>, const Cx *env) {
  return eval(L{}, env) * eval(R{}, env);
}

/// @brief Evaluate an expression type against a runtime environment.
///
/// This is the eager *tree* walk: a subexpression shared between branches is
/// recomputed at each occurrence. For maximal reuse on hot paths, evaluate the
/// expression as a DAG instead — see @ref numtracer::expr::eval_cse and @ref
/// numtracer::expr::eval_real.
/// @tparam E The expression type to evaluate.
/// @param env The runtime environment array; `env[Id]` supplies each `Var<Id>`.
/// @return The complex value of the expression.
template <class E> constexpr Cx evaluate(const Cx *env) {
  return eval(E{}, env);
}

} // namespace numtracer::expr
