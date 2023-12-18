use std::cmp::Ordering;

extern crate cairo;
extern crate rand;

fn draw_bounds(context: &cairo::Context, points: &[(f64, f64)], r: f64, g: f64, b: f64) {
    let mut min_x = points[0].0;
    let mut min_y = points[0].1;
    let mut max_x = points[0].0;
    let mut max_y = points[0].1;

    for (x, y) in points {
        if *x < min_x {
            min_x = *x
        };
        if *y < min_y {
            min_y = *y
        };
        if *x > max_x {
            max_x = *x
        };
        if *y > max_y {
            max_y = *y
        };
    }

    context.rectangle(min_x, min_y, max_x - min_x, max_y - min_y);
    context.set_source_rgb(r, g, b);
    context.set_line_width(0.01);
    context.stroke().unwrap();
}

fn sort_key(value: &(f64, f64), dimension: usize) -> f64 {
    if dimension % 2 == 0 {
        value.0
    } else {
        value.1
    }
}

fn main() {
    let surface = cairo::PdfSurface::new(100.0, 100.0, "out.pdf").unwrap();
    let context = cairo::Context::new(&surface).unwrap();

    use rand::Rng;
    let mut rng = rand::XorShiftRng::new_unseeded();

    let mut points = Vec::new();
    const NUM_POINTS: usize = 70001;
    const NUM_CHILD_NODES: usize = 25;

    println!("We want to generate a tree of {NUM_POINTS} elements with a maximum of {NUM_CHILD_NODES} children per node.");
    let tree_height = (NUM_POINTS as f64).log(NUM_CHILD_NODES as f64).ceil() as u32;
    println!("  This requires a tree height of {tree_height} ({NUM_CHILD_NODES}^{}={} < {} <= {NUM_CHILD_NODES}^{tree_height}={})",
        tree_height-1, NUM_CHILD_NODES.pow(tree_height-1), NUM_POINTS, NUM_CHILD_NODES.pow(tree_height));

    let max_num_children_per_node =
        (NUM_POINTS as f64).powf(1.0 / tree_height as f64).ceil() as usize;
    println!(
        "    With that height, we only need <={max_num_children_per_node} children per node (\
            {}^{tree_height}={} < {NUM_POINTS} <= {max_num_children_per_node}^{tree_height}={})",
        max_num_children_per_node - 1,
        (max_num_children_per_node - 1).pow(tree_height),
        max_num_children_per_node.pow(tree_height)
    );

    for _i in 0..NUM_POINTS {
        points.push((rng.next_f64() * 100.0, rng.next_f64() * 100.0));
    }

    for i in 0..=tree_height {
        println!(
            "On tree level {i}, we target having {} blocks",
            max_num_children_per_node.pow(i)
        );
    }

    let num_blocks = max_num_children_per_node.pow(tree_height - 1);
    let mut blocks = Vec::new();
    println!("For our first iteration, we target merging the {NUM_POINTS} points into exactly {num_blocks} blocks");
    let max_num_children_per_node = ((NUM_POINTS as f64) / (num_blocks as f64)).ceil() as usize;
    println!("This requires up to {max_num_children_per_node} nodes per block.");
    let mut start_index = 0;
    let mut remaining_points = NUM_POINTS;
    for block in 0..num_blocks {
        let remaining_blocks = num_blocks - block;
        //println!("Remaining: {remaining_points} points, {remaining_blocks} blocks");
        let num_points_in_block =
            match ((max_num_children_per_node - 1) * remaining_blocks).cmp(&remaining_points) {
                Ordering::Less => max_num_children_per_node,
                Ordering::Equal => max_num_children_per_node - 1,
                Ordering::Greater => panic!(),
            };

        println!("Block {block} has {num_points_in_block} points.");
        blocks.push((start_index, num_points_in_block));
        start_index += num_points_in_block;
        remaining_points -= num_points_in_block;
    }
    assert_eq!(remaining_points, 0);
    println!("Total {} blocks.", blocks.len());

    let num_slices = (num_blocks as f64).sqrt().ceil() as usize;
    let max_num_blocks_per_slice = ((num_blocks as f64) / (num_slices as f64)).ceil() as usize;
    println!("This requires partitioning the points into {num_slices} slices with up to {max_num_blocks_per_slice} blocks per slice");

    //let num_nodes_per_slice = NUM_POINTS as f64 / ((num_blocks as f64).sqrt());
    //let num_slices = (NUM_POINTS as f64 / num_nodes_per_slice).ceil() as usize ;
    //println!("#slices: {}", num_slices);

    //let num_nodes_per_slice = (NUM_POINTS as f64 / num_slices as f64).ceil() as usize;
    //println!("Nodes per slice: {}", num_slices *max_num_children_per_node);

    let mut start_block_index = 0;
    let mut remaining_blocks = num_blocks;

    points.sort_by(|lhs, rhs| sort_key(lhs, 0).partial_cmp(&sort_key(rhs, 0)).unwrap());

    for i in 0..num_slices {
        let remaining_slices = num_slices - i;
        let blocks_in_this_slice =
            match ((max_num_blocks_per_slice - 1) * remaining_slices).cmp(&remaining_blocks) {
                Ordering::Less => max_num_blocks_per_slice,
                Ordering::Equal => max_num_blocks_per_slice - 1,
                Ordering::Greater => panic!(),
            };

        let mut nodes_in_this_slice = 0;
        let this_slice_blocks =
            &blocks[start_block_index..start_block_index + blocks_in_this_slice];
        for (_, count) in this_slice_blocks {
            nodes_in_this_slice += count;
        }

        assert!(!this_slice_blocks.is_empty());
        let start_index = this_slice_blocks.first().unwrap().0;
        let (last_start, last_count) = this_slice_blocks.last().unwrap();
        let last_index = last_start + last_count;
        assert_eq!(last_index, start_index + nodes_in_this_slice);
        println!("Slice {i}, {blocks_in_this_slice} blocks, {nodes_in_this_slice} nodes");

        let slice = &mut points[start_index..start_index + nodes_in_this_slice];
        //draw_bounds(&context, &slice, 0.1, 0.1, (i as f64) *0.1);

        slice.sort_by(|lhs, rhs| sort_key(lhs, 1).partial_cmp(&sort_key(rhs, 1)).unwrap());

        for (block_start_index, nodes_in_this_block) in this_slice_blocks {
            let the_box = &points[*block_start_index..*block_start_index + *nodes_in_this_block];
            draw_bounds(&context, the_box, 0.1, (i as f64) * 0.1, 0.1);
        }
        start_block_index += blocks_in_this_slice;
        remaining_blocks -= blocks_in_this_slice;
    }

    let mut r: f64 = 0.0;
    for (x, y) in points {
        context.rectangle(x - 0.05, y - 0.05, 0.1, 0.1);
        context.set_line_width(0.1);
        context.set_source_rgb(r, 0.0, 0.0);
        r += 0.001;
        context.fill().unwrap();
    }
    println!("Created a total of {} blocks", start_block_index);
}
