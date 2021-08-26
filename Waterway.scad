union() {
    difference() {
        intersection() {
            import("temp/poly_1.stl");
            import("temp/waterway_1.stl");
            import("temp/terrain_low1.stl");
        };
        import("temp/road_3.stl");
        import("temp/trail_3.stl");
        union() {
            import("temp/waterbody_N1_1.stl");
            import("temp/waterbody_N1_2.stl");
            import("temp/waterbody_N1_3.stl");
            import("temp/waterbody_N1_4.stl");
            import("temp/waterbody_N1_5.stl"); 
            };
        };
    difference() {
        intersection() {
            import("temp/terrain_low1.stl");
            import("temp/poly_1.stl");
            import("temp/waterway_2.stl");
            };
        import("temp/road_3.stl");
        import("temp/trail_3.stl");
        import("temp/terrain_low2.stl");
        union() {
            import("temp/waterbody_N1_1.stl");
            import("temp/waterbody_N1_2.stl");
            import("temp/waterbody_N1_3.stl");
            import("temp/waterbody_N1_4.stl");
            import("temp/waterbody_N1_5.stl");
            };
        };
    };


