using Timekeepers, GLMakie

println("[warmup] Loading Timekeepers and GLMakie...")
println("[warmup] Building the TKApp window (compiles UI code)...")
app = TKApp(; size = (900, 620))

println("[warmup] Rendering once so GL/plotting paths get traced into the image...")
display(app)
sleep(2)

println("[warmup] Closing window.")
GLMakie.closeall()

println("[warmup] Done. These compiled paths will be baked into the sysimage.")
