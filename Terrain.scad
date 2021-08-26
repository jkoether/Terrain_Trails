difference() {
    intersection() {
        import("temp/terrain.stl");
        import("temp/poly_3.stl");
    };
    intersection() {
        import("temp/poly_2.stl");
        import("temp/road_3.stl");
    };
    intersection() {
        import("temp/poly_2.stl");
        import("temp/trail_3.stl");
    };
    intersection() {
        import("temp/poly_2.stl");
        difference() {
            import("temp/waterway_3.stl");
            union() {
                import("temp/waterbody_N1_1.stl");
                import("temp/waterbody_N1_2.stl");
                import("temp/waterbody_N1_3.stl");
                import("temp/waterbody_N1_4.stl");
                import("temp/waterbody_N1_5.stl"); 
            };
        };
    };
    intersection() {
        import("temp/poly_2.stl");
        union() {
            import("temp/waterbody_N2_1.stl");
            import("temp/waterbody_N2_2.stl");
            import("temp/waterbody_N2_3.stl");
            import("temp/waterbody_N2_4.stl");
            import("temp/waterbody_N2_5.stl"); 
            
        };
    };
};