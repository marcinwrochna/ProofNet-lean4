# ProofNet-lean4
This is a Lean 4.19.0 version of the ProofNet dataset, originally https://github.com/zhangir-azerbayev/ProofNet.

The original dataset (Lean 3) was released and maintained by Zhangir Azerbayev, Bartosz Piotrowski,  Hailey Schoelkopf, Edward W. Ayers, Dragomir Radev, and Jeremy Avigad (*ProofNet: Autoformalizing and Formally Proving Undergraduate-Level Mathematics*, [arxiv:2302.12433](https://arxiv.org/abs/2302.12433)), under an MIT license.

This repo is based on the port to Lean 4.6 by [Rahul Vishwakarma](https://github.com/rahul3613/ProofNet-lean4)'s (now updated to 4.18.0) and the [technical refactoring](https://github.com/njuyxw/Proofnet-lean4/compare/60efffb605ee07bf723db4fb8058129a7c8a89bb...main) by njuyxw (Xiao-Wen Yang 杨骁文), both also MIT-licensed.

This repo also fixes the formalization of putnam_exercise_2001_a5 (`∃! a : ℕ, ∃! n : ℕ` vs `∃! p : ℕ × ℕ`).
(Now also fixed in Rahul's version).
