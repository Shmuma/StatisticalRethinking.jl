ll = reshape([1,1,1,2], 2, 2)
r = waic(ll)
@test r.WAIC ≈ -4.740229013916554
@test r.lppd ≈ 2.620114506958277
@test r.penalty ≈ 0.25
@test r.std_err ≈ 0.5234209553714287

r = waic(ll, pointwise=true)
@test r.WAIC ≈ [-2.0, -2.740229013916555]
@test r.lppd ≈ [1.0, 1.6201145069582774]
@test r.penalty ≈ [0.0, 0.25]
@test r.std_err ≈ 0.5234209553714287
