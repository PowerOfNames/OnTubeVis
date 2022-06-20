#pragma once

// C++ STL
#include <vector>
#include <set>

// CGV framework core
#include <cgv/base/node.h>
#include <cgv/gui/event_handler.h>
#include <cgv/gui/provider.h>
#include <cgv/render/drawable.h>
#include <cgv/utils/stopwatch.h>

// CGV OpenGL lib
#include <cgv_gl/volume_renderer.h>

// CGV framework graphics utility
#include <cgv_glutil/application_plugin.h>
#include <cgv_glutil/box_wire_render_data.h>
#include <cgv_glutil/cone_render_data.h>
#include <cgv_glutil/color_map_editor.h>
#include <cgv_glutil/frame_buffer_container.h>
#include <cgv_glutil/navigator.h>
#include <cgv_glutil/radix_sort_4way.h>
#include <cgv_glutil/shader_library.h>
#include <cgv_glutil/sphere_render_data.h>
//#include <cgv_glutil/transfer_function_editor.h>

// local includes
#include "traj_loader.h"
#include "arclen_helper.h"
#include "demo.h" // interactive testbed helper classes and data
#include "attrib_handle_manager.h"
#include "voxelizer.h"
#include "ambient_occlusion_style.h"
#include "glyph_layer_manager.h"
#include "color_map_manager.h"
#include "textured_spline_tube_renderer.h"
#include "color_map_viewer.h"
#include "optix_integration.h"



using namespace cgv::render;

////
// Plugin definition

/// baseline visualization plugin for arbitrary trajectory data as tubes using the framework tube renderers and
/// trajectory loading facilities
class tubes :
	public cgv::base::argument_handler, // derive from argument handler to be able to process custom arguments
	public cgv::glutil::application_plugin	// derive from application plugin, which is a node, drawable, gui provider and event handler and can handle overlays
{
public:
	/// data layout for per-node attributes within the attribute render SSBO
	struct node_attribs {
		vec4 pos_rad;
		vec4 color;
		vec4 tangent;
	};

	cgv::type::DummyEnum voxel_grid_resolution;

protected:

	// ###############################
	// ### BEGIN: OptiX integration
	// ###############################

	// state
	struct
	{
		// use OptiX instead of OpenGL rasterization
		bool enabled = false;

		// subdivide into quadratic beziers for higher ray tracing performance
		bool subdivide = false;

		// OptiX device context
		OptixDeviceContext context = nullptr;

		// acceleration datastructure resources
		OptixTraversableHandle accelds = 0;
		CUdeviceptr            accelds_outbuf = 0;

		// curve tracing pipeline resources
		// - pipeline object
		OptixPipeline pipeline = nullptr;
		// - program groups
		OptixProgramGroup prg_hit = nullptr;
		OptixProgramGroup prg_miss = nullptr;
		OptixProgramGroup prg_raygen = nullptr;
		// - modules
		OptixModule mod_shading = nullptr;
		OptixModule mod_geom = nullptr;
		// - shader binding table
		OptixShaderBindingTable sbt;

		// device memory for our launch parameters
		CUdeviceptr params_buf = 0;

		// GL interop
		CUstream stream = nullptr;
		cuda_output_buffer<float4> outbuf_albedo;
		cuda_output_buffer<float1> outbuf_depth;
		texture tex_albedo, tex_depth;
	} optix;

	bool optix_update_accelds (void);
	bool optix_update_pipeline (void);

	void optix_destroy_accelds (void);
	void optix_destroy_pipeline (void);

	void optix_draw_trajectories (context &ctx);

	// ###############################
	// ###  END:  OptiX integration
	// ###############################


	cgv::glutil::color_map_editor_ptr cm_editor_ptr;
	cgv::glutil::color_map_editor_ptr tf_editor_ptr;
	cgv::data::ref_ptr<cgv::glutil::navigator> navigator_ptr;
	cgv::data::ref_ptr <color_map_viewer> cm_viewer_ptr;

	enum GridMode {
		GM_NONE = 0,
		GM_COLOR = 1,
		GM_NORMAL = 2,
		GM_COLOR_NORMAL = 3
	};

	struct grid_parameters {
		vec2 scaling;
		float thickness;
		float blend_factor;
	};

	GridMode grid_mode;
	rgba grid_color;
	cgv::type::DummyEnum grid_normal_settings;
	bool grid_normal_inwards;
	bool grid_normal_variant;
	float normal_mapping_scale;
	std::vector<grid_parameters> grids;
	bool enable_fuzzy_grid;

	/// global tube render settings
	struct {
		bool use_curvature_correction = true;
		float length_scale = 1.0;
		float antialias_radius = 0.5f;
	} general_settings;
	
	/// shader defines for the deferred shading pass
	shader_define_map tube_shading_defines;

	/// store a pointer to the view for fast access
	view* view_ptr = nullptr;
	
	/// store the current OpenGL viewport configuration
	GLint viewport[4];

	/// path of the dataset to load - can be either a directory or a single file
	std::string datapath;

	/// misc configurable fields
	struct {
		/// proxy for controlling fltk_gl_view::instant_redraw
		bool instant_redraw_proxy = false;

		/// proxy for controlling context::enable_vsynch through fltk_gl_view
		bool vsync_proxy = true;

		/// proxy for controlling stereo_view_interactor::fix_view_up_dir
		bool fix_view_up_dir_proxy = false;
	} misc_cfg;

	/// drag-n-drop state fields
	struct {
		/// current mouse position
		ivec2 pos;

		/// current drag-n-drop string
		std::string text;

		/// list of filenames extracted from @ref #text
		std::vector<std::string> filenames;
	} dnd;

	/// dataset state fields
	struct {
		/// set of filepaths for loading
		std::set<std::string> files;

		/// generated demo dataset
		std::vector<demo::trajectory> demo_trajs;
	} dataset;

	// TODO: comment and cleanup members
	cgv::glutil::frame_buffer_container fbc;
	cgv::glutil::shader_library shaders;
	volume_render_style vstyle;
	cgv::glutil::gl_color_map volume_tf;

	bool show_bbox = true;
	bool show_wireframe_bbox = true;
	cgv::render::box_render_style bbox_style;
	cgv::glutil::box_render_data<> bbox_rd;
	cgv::glutil::box_wire_render_data<> bbox_wire_rd;

	/// rendering state fields
	struct {
		/// render style for the textured spline tubes
		textured_spline_tube_render_style style;
		
		/// render data generated by the trajectory manager
		const traj_manager<float>::render_data *data;

		/// segment-wise arclength approximations (set of 4 cubic bezier curves returning global
		/// trajectory arclength at the segment, packed into the columns of a 4x4 matrix)
		arclen::parametrization arclen_data;

		/// GPU-side storage buffer mirroring the \ref #arclen_data .
		vertex_buffer arclen_sbo;

		/// GPU-side render attribute buffer.
		vertex_buffer render_sbo;

		/// GPU-side storage buffers storing independently sampled attribute data.
		std::vector<vertex_buffer> attribs_sbos;

		/// GPU-side storage buffers indexing the independently sampled attributes per tube segment.
		std::vector<vertex_buffer> aindex_sbos;

		/// shared attribute array manager used by both renderers
		attribute_array_manager aam;

		/// the gpu sorter used to reorder the indices according to their corresponding segment visibility order
		cgv::glutil::gpu_sorter* sorter = nullptr;
	} render;

	/// trajectory manager
	traj_manager<float> traj_mgr;

	/// attribute handle manager
	attrib_handle_manager<float> ah_mgr;

	/// glyph layer manager
	glyph_layer_manager glyph_layer_mgr;

	/// color map manager
	color_map_manager color_map_mgr;

	/// benchmark state fields
	struct {
		/// whether a benchmark run is requested
		bool requested = false;
		/// whether a benchmark is currently running
		bool running = false;
		/// timer to count elapsed time
		cgv::utils::stopwatch timer;
		/// counter for rendered frames
		unsigned total_frames;
		/// store last seconds since the start of the run
		double last_seconds_since_start;
	} benchmark;

	/// the different debug render modes
	enum DebugRenderMode {
		DRM_NONE,
		DRM_NODES,
		DRM_SEGMENTS,
		DRM_NODES_SEGMENTS
	};

	/// debug state fields
	struct {
		DebugRenderMode render_mode = DRM_NONE;

		/// debug render data and styles
		cgv::glutil::sphere_render_data<> node_rd;
		cgv::glutil::cone_render_data<> segment_rd;
		sphere_render_style node_rs;
		cone_render_style segment_rs;

		/// whether to higlight individual segments in the textured spline tube renderer
		bool highlight_segments = false;

		/// whether to sort the segments, which is used to boost performance together with conservative depth testing
		bool sort = true;
		/// whether to foirce the initial draw order of segments as defined in the data set (overrides sort setting)
		bool force_initial_order = false;
		/// whether to limit the render count
		bool limit_render_count = false;
		/// percentage of rendered segments
		float render_percentage = 1.0f;
		/// amount of rendered segments
		size_t render_count = 0;
		/// total segment count
		size_t segment_count = 0;

		double far_extent_factor = 0.8;
		double near_extent_factor = 0.3;
		bool near_view = false;
	} debug;

	bool benchmark_mode = false;
	bool benchmark_mode_setup = false;
	
	/// members for rendering eye position and direction used to test sorting
	cgv::glutil::sphere_render_data<> srd;
	vec3 test_eye = vec3(5.0f, 0.5f, 5.0f);
	vec3 test_dir = vec3(0.0f, 0.0f, -1.0f);

	/// file handling fields
	struct {
		/// file name of loaded layer configuration
		std::string file_name;
		/// file name of to be saved layer configuration (used to trigger save action)
		std::string save_file_name;
		/// track whether the current configuration has unsaved changes
		bool has_unsaved_changes = false;
	} fh;

	bool voxelize_gpu = false;
	bool show_volume = false;
	box3 bbox;
	texture density_tex;
	texture tf_tex;

	voxelizer density_volume;
	ambient_occlusion_style ao_style, ao_style_bak; // the latter is used to restore defaults after demo data is unloaded

	glyph_layer_manager::configuration glyph_layers_config;
	bool include_hidden_glyphs = false;
	unsigned max_glyph_count = 10;
	
	void reload_shader();
	bool save_layer_configuration(const std::string& file_name);
	bool read_layer_configuration(const std::string& file_name);

	void update_glyph_layer_manager(void);
	void glyphs_out_of_date(bool state);
	bool compile_glyph_attribs(void);
	double change_time = 0.0;
	double recalc_delay = 0.2;
	bool has_changed = false;
	void timer_event(double, double);

	void set_view(void);
	void update_grid_ratios(void);
	void update_attribute_bindings(void);
	void update_debug_attribute_bindings(void);
	void calculate_bounding_box(void);

	void create_density_volume(context& ctx, unsigned resolution);

	/// draw methods
	void draw_dnd(context& ctx);
	void draw_trajectories(context& ctx);
	void draw_density_volume(context& ctx);

	/// helper methods
	shader_define_map build_tube_shading_defines();

	void create_vec3_gui(const std::string& name, vec3& value, float min = 0.0f, float max = 1.0f);

public:
	tubes();
	~tubes();

	std::string get_type_name() const { return "tubes"; }
	void handle_args(std::vector<std::string> &args);

	void clear(context& ctx);

	bool self_reflect(cgv::reflect::reflection_handler& rh);
	void stream_help(std::ostream& os);
	void stream_stats(std::ostream& os) {}

	bool handle_event(cgv::gui::event& e);
	void on_set(void* member_ptr);
	bool on_exit_request();

	bool init(context& ctx);
	void init_frame(context& ctx);
	void draw(context& ctx);

	void create_gui();
};
