#include <iostream>
#include <vector>

#include <sycl/sycl.hpp>

int main(int argc, char** argv){
	if(argc > 1) {
		std::cerr << "Invalid argument: \"" << argv[1] << "\"" << std::endl;

		return EXIT_FAILURE;
	}

	std::cout << "Hello, world!" << std::endl;

	sycl::queue q;
	
	std::size_t input_size = 1024;

	std::vector<int> input(input_size);

	for(int i = 0; i < input.size(); i++) input[i] = i;

	sycl::buffer<int> buff{input.data(), sycl::range<1>{input_size}};

	constexpr size_t Group_size = 128;

	q.submit([&](sycl::handler& cgh){
		auto data_accessor = buff.get_access<sycl::access::mode::read_write>(cgh);
		cgh.parallel<class Kernel>(sycl::range<1>{input_size / Group_size}, sycl::range<1>{Group_size}, 
		[=](auto grp){
			// Outside of distribute_items(), the degree of parallelism is implementation-defined.
			// the implementation can use whatever is most efficient for hardware/backend.
			// In hipSYCL CPU, this would be executed by a single thread on CPU
			// and Group_size threads on hipSYCL GPU
			// Information about the position in the physical iteration space can be obtained
			// using grp.get_physical_local_id() and grp.get_physical_local_range().

			// sycl::memory_environment() can be used to allocate local memory 
			// (of compile-time size) as well as private memory that is persistent across
			// multiple distribute_items() calls.
			// Of course, local accessors can also be used.
			sycl::memory_environment(grp, 
				sycl::require_local_mem<int[Group_size]>(),
				// the requested private memory is not used in this example,
				// and only here to showcase how to request private memory.
				sycl::require_private_mem<int>(),
				[&](auto& scratch, auto& private_mem){
				// the arguments passed to the lambda function corresponds to the require arguments before.
				// private_mem is allocated in private memory of the _logical_ work item
				// and of type sycl::s_private_memory<T, decltype(grp)>&.
				// scratch is a reference to the requested int[Group_size] array.

				// Variables not explicitly requested as local or private memory 
				// will be allocated in private memory of the _physical_ work item
				// (see the for loop below)

				// `distribute_items` distributes the logical, user-provided iteration space across the physical one. 
				sycl::distribute_items(grp, [&](sycl::s_item<1> idx){
						scratch[idx.get_local_id(grp, 0)] = data_accessor[idx.get_global_id(0)]; 
				});
				// Instead of an explicit group_barrier, we could also use the
				// blocking distribute_items_and_wait()
				sycl::group_barrier(grp);

				// Can execute code e.g. for a single item of a subgroup:
				sycl::distribute_groups(grp, [&](auto subgroup){
					sycl::single_item(subgroup, [&](){
						// ...
					});
				});

				// Variables inside the parallel scope that are not explicitly local or private memory
				// are allowed, if they are not modified from inside `distribute_items()` scope.
				// The SYCL implementation will allocate those in private memory of the physical item,
				// so they will always be efficient. This implies that the user should not attempt to assign values
				// per logical work item, since they are allocated per physical item.
				for(int i = Group_size / 2; i > 0; i /= 2){
					// The *_and_wait variants of distribute_groups and distribute_items
					// invoke a group_barrier at the end.
					sycl::distribute_items_and_wait(grp, 
						[&](sycl::s_item<1> idx){
						size_t lid = idx.get_innermost_local_id(0);
						if(lid < i)
							scratch[lid] += scratch[lid+i];
					});
				}
				
				sycl::single_item(grp, [&](){
					data_accessor[grp.get_group_id(0)*Group_size] = scratch[0];
				});
			});
		});
	});

	// Verify results on host
	auto host_acc = buff.get_access<sycl::access::mode::read>();

	for(int grp = 0; grp < input_size/Group_size; grp++){
		int host_result = 0;

		for(int i = grp * Group_size; i < (grp + 1) * Group_size; i++) host_result += i;

		if(host_result != host_acc[grp*Group_size]) std::cout << "Wrong result, got " << host_acc[grp * Group_size] << ", expected " << host_result << std::endl;
	}

	return EXIT_SUCCESS;
}
