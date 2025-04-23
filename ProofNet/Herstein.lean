import Mathlib

open Fintype Set Real Ideal Polynomial
open scoped BigOperators
noncomputable section

theorem herstein_exercise_2_1_18 {G : Type*} [Group G]
  [Fintype G] (hG2 : Even (card G)) :
  ∃ (a : G), a ≠ 1 ∧ a = a⁻¹ :=
sorry

theorem herstein_exercise_2_1_21 (G : Type*) [Group G] [Fintype G]
  (hG : card G = 5) :
  Nonempty (CommGroup G) :=
sorry

theorem herstein_exercise_2_1_26 {G : Type*} [Group G]
  [Fintype G] (a : G) : ∃ (n : ℕ), a ^ n = 1 :=
sorry

theorem herstein_exercise_2_1_27 {G : Type*} [Group G]
  [Fintype G] : ∃ (m : ℕ), ∀ (a : G), a ^ m = 1 :=
sorry

theorem herstein_exercise_2_2_3 {G : Type*} [Group G]
  {P : ℕ → Prop} {hP : P = λ i => ∀ a b : G, (a*b)^i = a^i * b^i}
  (hP1 : ∃ n : ℕ, P n ∧ P (n+1) ∧ P (n+2)) : Nonempty (CommGroup G) :=
sorry

theorem herstein_exercise_2_2_5 {G : Type*} [Group G]
  (h : ∀ (a b : G), (a * b) ^ 3 = a ^ 3 * b ^ 3 ∧ (a * b) ^ 5 = a ^ 5 * b ^ 5) :
  Nonempty (CommGroup G) :=
sorry

theorem herstein_exercise_2_2_6c {G : Type*} [Group G] {n : ℕ} (hn : n > 1)
  (h : ∀ (a b : G), (a * b) ^ n = a ^ n * b ^ n) :
  ∀ (a b : G), (a * b * a⁻¹ * b⁻¹) ^ (n * (n - 1)) = 1 :=
sorry

theorem herstein_exercise_2_3_17 {G : Type*} [Mul G] [Group G] (a x : G) :
  centralizer {x⁻¹*a*x} =
  (λ g : G => x⁻¹*g*x) '' (centralizer {a}) :=
sorry

theorem herstein_exercise_2_3_16 {G : Type*} [Group G]
  (hG : ∀ H : Subgroup G, H = ⊤ ∨ H = ⊥) :
  IsCyclic G ∧ ∃ (p : ℕ) (Fin : Fintype G), Nat.Prime p ∧ @card G Fin = p :=
sorry

theorem herstein_exercise_2_4_36 {a n : ℕ} (h : a > 1) :
  n ∣ (a ^ n - 1).totient :=
sorry

theorem herstein_exercise_2_5_23 {G : Type*} [Group G]
  (hG : ∀ (H : Subgroup G), H.Normal) (a b : G) :
  ∃ (j : ℤ) , b*a = a^j * b:=
sorry

theorem herstein_exercise_2_5_30 {G : Type*} [Group G] [Fintype G]
  {p m : ℕ} (hp : Nat.Prime p) (hp1 : ¬ p ∣ m) (hG : card G = p*m)
  {H : Subgroup G} [Fintype H] [H.Normal] (hH : card H = p):
  Subgroup.Characteristic H :=
sorry

theorem herstein_exercise_2_5_31 {G : Type*} [CommGroup G] [Fintype G]
  {p m n : ℕ} (hp : Nat.Prime p) (hp1 : ¬ p ∣ m) (hG : card G = p^n*m)
  {H : Subgroup G} [Fintype H] (hH : card H = p^n) :
  Subgroup.Characteristic H :=
sorry

theorem herstein_exercise_2_5_37 (G : Type*) [Group G] [Fintype G]
  (hG : card G = 6) (hG' : IsEmpty (CommGroup G)) :
  Nonempty (G ≃* Equiv.Perm (Fin 3)) :=
sorry

theorem herstein_exercise_2_5_43 (G : Type*) [Group G] [Fintype G]
  (hG : card G = 9) :
  Nonempty (CommGroup G) :=
sorry

theorem herstein_exercise_2_5_44 {G : Type*} [Group G] [Fintype G] {p : ℕ}
  (hp : Nat.Prime p) (hG : card G = p^2) :
  ∃ (N : Subgroup G) (Fin : Fintype N), @card N Fin = p ∧ N.Normal :=
sorry

theorem herstein_exercise_2_5_52 {G : Type*} [Group G] [Fintype G]
  (φ : G ≃* G) {I : Finset G} (hI : ∀ x ∈ I, φ x = x⁻¹)
  (hI1 : (0.75 : ℚ) * card G ≤ card I) :
  ∀ x : G, φ x = x⁻¹ ∧ ∀ x y : G, x*y = y*x :=
sorry

theorem herstein_exercise_2_6_15 {G : Type*} [CommGroup G] {m n : ℕ}
  (hm : ∃ (g : G), orderOf g = m)
  (hn : ∃ (g : G), orderOf g = n)
  (hmn : m.Coprime n) :
  ∃ (g : G), orderOf g = m * n :=
sorry

theorem herstein_exercise_2_7_7 {G : Type*} [Group G] {G' : Type*} [Group G']
  (φ : G →* G') (N : Subgroup G) [N.Normal] :
  (Subgroup.map φ N).Normal  :=
sorry

theorem herstein_exercise_2_8_12 {G H : Type*} [Fintype G] [Fintype H]
  [Group G] [Group H] (hG : card G = 21) (hH : card H = 21)
  (hG1 : IsEmpty (CommGroup G)) (hH1 : IsEmpty (CommGroup H)) :
  Nonempty (G ≃* H) :=
sorry

theorem herstein_exercise_2_8_15 {G H: Type*} [Fintype G] [Group G] [Fintype H]
  [Group H] {p q : ℕ} (hp : Nat.Prime p) (hq : Nat.Prime q)
  (h : p > q) (h1 : q ∣ p - 1) (hG : card G = p*q) (hH : card G = p*q) :
  Nonempty (G ≃* H) :=
sorry

theorem herstein_exercise_2_9_2 {G H : Type*} [Fintype G] [Fintype H] [Group G]
  [Group H] (hG : IsCyclic G) (hH : IsCyclic H) :
  IsCyclic (G × H) ↔ (card G).Coprime (card H) :=
sorry

theorem herstein_exercise_2_10_1 {G : Type*} [Group G] (A : Subgroup G)
  [A.Normal] {b : G} (hp : Nat.Prime (orderOf b)) :
  A ⊓ (Subgroup.closure {b}) = ⊥ :=
sorry

theorem herstein_exercise_2_11_6 {G : Type*} [Group G] {p : ℕ} (hp : Nat.Prime p)
  {P : Sylow p G} (hP : P.Normal) :
  ∀ (Q : Sylow p G), P = Q :=
sorry

theorem herstein_exercise_2_11_7 {G : Type*} [Group G] {p : ℕ} (hp : Nat.Prime p)
  {P : Sylow p G} (hP : P.Normal) :
  Subgroup.Characteristic (P : Subgroup G) :=
sorry

theorem herstein_exercise_2_11_22 {p : ℕ} {n : ℕ} {G : Type*} [Fintype G]
  [Group G] (hp : Nat.Prime p) (hG : card G = p ^ n) {K : Subgroup G}
  [Fintype K] (hK : card K = p ^ (n-1)) :
  K.Normal :=
sorry

theorem herstein_exercise_3_2_21 {α : Type*} [Fintype α] {σ τ: Equiv.Perm α}
  (h1 : ∀ a : α, σ a = a ↔ τ a ≠ a) (h2 : τ ∘ σ = id) :
  σ = 1 ∧ τ = 1 :=
sorry

theorem herstein_exercise_4_1_19 : Infinite {x : Quaternion ℝ | x^2 = -1} :=
sorry

theorem herstein_exercise_4_1_34 : Nonempty $ Equiv.Perm (Fin 3) ≃* Matrix.GeneralLinearGroup (Fin 2) (ZMod 2) :=
sorry

theorem herstein_exercise_4_2_5 {R : Type*} [Ring R]
  (h : ∀ x : R, x ^ 3 = x) : Nonempty (CommRing R) :=
sorry

theorem herstein_exercise_4_2_6 {R : Type*} [Ring R] (a x : R)
  (h : a ^ 2 = 0) : a * (a * x + x * a) = (x + x * a) * a :=
sorry

theorem herstein_exercise_4_2_9 {p : ℕ} (hp : Nat.Prime p) (hp1 : Odd p) :
  ∃ (a b : ℤ), (a / b : ℚ) = ∑ i ∈ Finset.range p, 1 / (i + 1) → ↑p ∣ a :=
sorry

theorem herstein_exercise_4_3_1 {R : Type*} [CommRing R] (a : R) :
  ∃ I : Ideal R, {x : R | x*a=0} = I :=
sorry

theorem herstein_exercise_4_3_25 (I : Ideal (Matrix (Fin 2) (Fin 2) ℝ)) :
  I = ⊥ ∨ I = ⊤ :=
sorry

theorem herstein_exercise_4_4_9 (p : ℕ) (hp : Nat.Prime p) :
  (∃ S : Finset (ZMod p), S.card = (p-1)/2 ∧ ∃ x : ZMod p, x^2 = p) ∧
  (∃ S : Finset (ZMod p), S.card = (p-1)/2 ∧ ¬ ∃ x : ZMod p, x^2 = p) :=
sorry

theorem herstein_exercise_4_5_16 {p n: ℕ} (hp : Nat.Prime p)
  {q : Polynomial (ZMod p)} (hq : Irreducible q) (hn : q.degree = n) :
  ∃ is_fin : Fintype $ Polynomial (ZMod p) ⧸ span ({q} : Set (Polynomial $ ZMod p)),
  @card (Polynomial (ZMod p) ⧸ span {q}) is_fin = p ^ n ∧
  IsField (Polynomial $ ZMod p):=
sorry

theorem herstein_exercise_4_5_23 {p q: Polynomial (ZMod 7)}
  (hp : p = X^3 - 2) (hq : q = X^3 + 2) :
  Irreducible p ∧ Irreducible q ∧
  (Nonempty $ Polynomial (ZMod 7) ⧸ span ({p} : Set $ Polynomial $ ZMod 7) ≃+*
  Polynomial (ZMod 7) ⧸ span ({q} : Set $ Polynomial $ ZMod 7)) :=
sorry

theorem herstein_exercise_4_5_25 {p : ℕ} (hp : Nat.Prime p) :
  Irreducible (∑ i : Finset.range p, X ^ p : Polynomial ℚ) :=
sorry

theorem herstein_exercise_4_6_2 : Irreducible (X^3 + 3*X + 2 : Polynomial ℚ) :=
sorry

theorem herstein_exercise_4_6_3 :
  Infinite {a : ℤ | Irreducible (X^7 + 15*X^2 - 30*X + (a : Polynomial ℚ) : Polynomial ℚ)} :=
sorry

theorem herstein_exercise_5_1_8 {p m n: ℕ} {F : Type*} [Field F]
  (hp : Nat.Prime p) (hF : CharP F p) (a b : F) (hm : m = p ^ n) :
  (a + b) ^ m = a^m + b^m :=
sorry

theorem herstein_exercise_5_2_20 {F V ι: Type*} [Infinite F] [Field F]
  [AddCommGroup V] [Module F V] {u : ι → Submodule F V}
  (hu : ∀ i : ι, u i ≠ ⊤) :
  (⋃ i : ι, (u i : Set V)) ≠ ⊤ :=
sorry

theorem herstein_exercise_5_3_7 {K : Type*} [Field K] {F : Subfield K}
  {a : K} (ha : IsAlgebraic F (a ^ 2)) : IsAlgebraic F a :=
sorry

theorem herstein_exercise_5_3_10 : IsAlgebraic ℚ (cos (Real.pi / 180)) :=
sorry

theorem herstein_exercise_5_4_3 {a : ℂ} {p : ℂ → ℂ}
  (hp : p = λ (x : ℂ) => x^5 + sqrt 2 * x^3 + sqrt 5 * x^2 + sqrt 7 * x + 11)
  (ha : p a = 0) :
  ∃ p : Polynomial ℂ , p.degree < 80 ∧ a ∈ p.roots ∧
  ∀ n : p.support, ∃ a b : ℤ, p.coeff n = a / b :=
sorry

theorem herstein_exercise_5_5_2 : Irreducible (X^3 - 3*X - 1 : Polynomial ℚ) :=
sorry

theorem herstein_exercise_5_6_14 {p m n: ℕ} (hp : Nat.Prime p) {F : Type*}
  [Field F] [CharP F p] (hm : m = p ^ n) :
  card (rootSet (X ^ m - X : Polynomial F) F) = m :=
sorry
