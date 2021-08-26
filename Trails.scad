
union() {
    difference() {
        intersection() {
            import("temp/poly_1.stl");
            import("temp/trail_1.stl");
            import("temp/terrain.stl");
        };
    import("temp/road_3.stl");
    };

    difference() {
        intersection() {
            import("temp/terrain.stl");
            import("temp/poly_1.stl");
            import("temp/trail_2.stl");
            };
        import("temp/terrain_low2.stl");
        import("temp/road_3.stl");
    };
};

