beta_values = 10 .^ range(log10(0.1), log10(10), length=101)

# Define tau values
tau_values = [1.0, 0.1, 0.01, 0.001, 0.0001]

# Convert arrays to comma-separated strings
tau_str = join(tau_values, ",")
beta_str = join(beta_values, ",")

# Write to config file
config_content = """# Tau values (comma-separated)
$tau_str
# Beta values (comma-separated)
$beta_str"""

# Write to file
open("tau_beta_config.txt", "w") do file
    write(file, config_content)
end

println("Config file created with:")
println("Tau values: $tau_str")
println("Beta values: $(length(beta_values)) values from $(beta_values[1]) to $(beta_values[end])")
println("File saved as tau_beta_config.txt")