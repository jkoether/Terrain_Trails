//import("temp/terrain.stl");
union() {
    intersection() {
            import("temp/poly_1.stl");
            import("temp/road_1.stl");
            import("temp/terrain.stl");
        };

    difference() {
        intersection() {
            import("temp/terrain.stl");
            import("temp/poly_1.stl");
            import("temp/road_2.stl");
            };

        import("temp/terrain_low2.stl");
    };
};

