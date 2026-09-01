import GalerkinToolkit as GT

function create_cartesian_mesh(domain, num_elements_per_dir)
    mesh = GT.cartesian_mesh(domain,num_elements_per_dir)
    return mesh
end


