namespace Garcia_vis_helper {

std::string data_root = "/scratch4/data/Garcia/";
std::string mesh_base = data_root + "microstructure/";
std::string stat_base = data_root + "Stat/";
std::string stress_base = data_root + "Stress/";
std::string movie_base = data_root + "Movies/";

dataset_info textured_info, untextured_info, mc_info_00, mc_info_06, mc_info_10_09, mc_info_10_05;

void set_paths()
{
    textured_info.mesh_base = mesh_base + "TexBNKT_190_190_37";
    textured_info.stress = stress_base + "vtk/trace-Estress_TexturedBNKT_b09.vtk";
    textured_info.stress_nrrd = stress_base + "nrrd/trace-Estress_TexturedBNKT_b09.nrrd";
    textured_info.dfield = data_root + "TexBNKT_revEngBNKT_b09_0/Dfield-scaled.nrrd";
    textured_info.id_to_tags = stat_base + "Textured_b09-ids_to_tags.nrrd";
    textured_info.stress_span = stat_base + "Textured_b09-grain-span.nrrd";
    textured_info.movie_path = "Textured_b09/stress/";
    textured_info.dfield = data_root + "TexBNKT_revEngBNKT_b09_0/Dfield-scaled.nrrd";
    textured_info.position = nvis::vec3(-38.2382, -20.2557, 27.1584);
    textured_info.focal_point = nvis::vec3(17.1675, 19.0469, 3.32833);
    textured_info.up = nvis::vec3(0.254529, 0.213122, 0.943289);
    textured_info.near = 17.6447;
    textured_info.far = 140.982;
    textured_info.step = nvis::vec3(0.2, 0.2, 0.4);
    textured_info.valid_bounds.min() = nvis::vec3(10., 10., 5.);
    textured_info.valid_bounds.max() = nvis::vec3(179., 179., 31.);
    
    untextured_info.mesh_base = mesh_base + "UnTexBNKT_01_140_140_70";
    untextured_info.stress = stress_base + "vtk/trace-Estress_UntexturedBNKT_b09.vtk";
    untextured_info.stress_nrrd = stress_base + "nrrd/trace-Estress_UntexturedBNKT_b09.nrrd";
    untextured_info.dfield = data_root + "UnTexBNKT_revEngBNKT_3m_b09_0/Dfield-scaled.nrrd";
    untextured_info.id_to_tags = stat_base + "Untextured_b09-ids_to_tags.nrrd";
    untextured_info.stress_span = stat_base + "Untextured_b09-grain-span.nrrd";
    untextured_info.movie_path = "Untextured_b09/stress/";
    untextured_info.position = nvis::vec3(-11.7274, -4.0038, 11.1954);
    untextured_info.focal_point = nvis::vec3(15.4841, 11.0543, -3.40036);
    untextured_info.up = nvis::vec3(0.352451, 0.239877, 0.904565);
    untextured_info.near = 5.28852;
    untextured_info.far = 39.2914;
    untextured_info.step = nvis::vec3(0.07, 0.07, 0.1);
    untextured_info.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    untextured_info.valid_bounds.max() = nvis::vec3(129., 129., 59.);
    
    mc_info_00.mesh_base = mesh_base + "MC_588grains_110cubed";
    mc_info_00.stress = stress_base + "vtk/trace-Estress_CG_588Grains_r00b09.vtk";
    mc_info_00.stress_nrrd = stress_base + "nrrd/trace-Estress_CG_588Grains_r00b09.nrrd";
    mc_info_00.dfield = data_root + "MC588GrainsRevEngBNKT_r00b09_0/Dfield-scaled.nrrd";
    mc_info_00.id_to_tags = stat_base + "CG_588Grains_r00b09-ids_to_tags.nrrd";
    mc_info_00.stress_span = stat_base + "CG_588Grains_r00b09-grain-span.nrrd";
    mc_info_00.movie_path = "MC_588grains/r00b09/stress/";
    mc_info_00.position = nvis::vec3(-324.175, -118.824, 264.273);
    mc_info_00.focal_point = nvis::vec3(38.2538, 78.0644, 115.621);
    mc_info_00.up = nvis::vec3(0.286275, 0.183103, 0.940489);
    mc_info_00.near = 157.65;
    mc_info_00.far = 961.655;
    mc_info_00.step = nvis::vec3(2, 2, 2);
    mc_info_00.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    mc_info_00.valid_bounds.max() = nvis::vec3(99., 99., 99.);
    
    mc_info_06.mesh_base = mesh_base + "MC_588grains_110cubed";
    mc_info_06.stress = stress_base + "vtk/trace-Estress_CG_588Grains_r06b09.vtk";
    mc_info_06.stress_nrrd = stress_base + "nrrd/trace-Estress_CG_588Grains_r06b09.nrrd";
    mc_info_06.dfield = data_root + "MC588GrainsRevEngBNKT_r06b09_0/Dfield-scaled.nrrd";
    mc_info_06.id_to_tags = stat_base + "CG_588Grains_r06b09-ids_to_tags.nrrd";
    mc_info_06.stress_span = stat_base + "CG_588Grains_r06b09-grain-span.nrrd";
    mc_info_06.movie_path = "MC_588grains/r06b09/stress/";
    mc_info_06.position = nvis::vec3(-324.175, -118.824, 264.273);
    mc_info_06.focal_point = nvis::vec3(38.2538, 78.0644, 115.621);
    mc_info_06.up = nvis::vec3(0.286275, 0.183103, 0.940489);
    mc_info_06.near = 157.65;
    mc_info_06.far = 961.655;
    mc_info_06.step = nvis::vec3(2, 2, 2);
    mc_info_06.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    mc_info_06.valid_bounds.max() = nvis::vec3(99., 99., 99.);
    
    mc_info_10_09.mesh_base = mesh_base + "MC_588grains_110cubed";
    mc_info_10_09.stress = stress_base + "vtk/trace-Estress_CG_588Grains_r10b09.vtk";
    mc_info_10_09.stress_nrrd = stress_base + "nrrd/trace-Estress_CG_588Grains_r10b09.nrrd";
    mc_info_10_09.dfield = data_root + "MC588GrainsRevEngBNKT_r10b09_0/Dfield-scaled.nrrd";
    mc_info_10_09.id_to_tags = stat_base + "CG_588Grains_r10b09-ids_to_tags.nrrd";
    mc_info_10_09.stress_span = stat_base + "CG_588Grains_r10b09-grain-span.nrrd";
    mc_info_10_09.movie_path = "MC_588grains/r10b09/stress/";
    mc_info_10_09.position = nvis::vec3(-324.175, -118.824, 264.273);
    mc_info_10_09.focal_point = nvis::vec3(38.2538, 78.0644, 115.621);
    mc_info_10_09.up = nvis::vec3(0.286275, 0.183103, 0.940489);
    mc_info_10_09.near = 157.65;
    mc_info_10_09.far = 961.655;
    mc_info_10_09.step = nvis::vec3(2, 2, 2);
    mc_info_10_09.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    mc_info_10_09.valid_bounds.max() = nvis::vec3(99., 99., 99.);
    
    mc_info_10_05.mesh_base = mesh_base + "MC_588grains_110cubed";
    mc_info_10_05.stress = stress_base + "vtk/trace-Estress_CG_588Grains_r10b05.vtk";
    mc_info_10_05.stress_nrrd = stress_base + "nrrd/trace-Estress_CG_588Grains_r10b05.nrrd";
    mc_info_10_05.dfield = data_root + "MC588GrainsRevEngBNKT_r10b05_0/Dfield-scaled.nrrd";
    mc_info_10_05.id_to_tags = stat_base + "CG_588Grains_r10b05-ids_to_tags.nrrd";
    mc_info_10_05.stress_span = stat_base + "CG_588Grains_r10b05-grain-span.nrrd";
    mc_info_10_05.movie_path = "MC_588grains/r10b05/stress/";
    mc_info_10_05.position = nvis::vec3(-324.175, -118.824, 264.273);
    mc_info_10_05.focal_point = nvis::vec3(38.2538, 78.0644, 115.621);
    mc_info_10_05.up = nvis::vec3(0.286275, 0.183103, 0.940489);
    mc_info_10_05.near = 157.65;
    mc_info_10_05.far = 961.655;
    mc_info_10_05.step = nvis::vec3(2, 2, 2);
    mc_info_10_05.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    mc_info_10_05.valid_bounds.max() = nvis::vec3(99., 99., 99.);
}

}


