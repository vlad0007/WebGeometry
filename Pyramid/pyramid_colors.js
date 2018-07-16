function facet_colors()
{
	var ind = 0;

	// table
	colors[ind] = new THREE.Color("rgb(200, 200, 200)"); // 0, 3, 2, 1, 0,    грань 0
	ind++;
	
	// crown
	colors[ind] = new THREE.Color("rgb(250, 190, 190)"); // 0, 4, 7, 3, 0,    грань 1
	ind++;
	colors[ind] = new THREE.Color("rgb(190, 250, 190)"); // 1, 5, 4, 0, 1,    грань 2
	ind++;
	colors[ind] = new THREE.Color("rgb(190, 190, 250)"); // 2, 6, 5, 1, 2,    грань 3
	ind++;
	colors[ind] = new THREE.Color("rgb(150, 250, 250)"); // 3, 7, 6, 2, 3,    грань 4
	ind++;
	
	// girdle
	colors[ind] = new THREE.Color("rgb(200, 200, 200)"); // 4, 8, 11, 7, 4,   грань 5
	ind++;
	colors[ind] = new THREE.Color("rgb(170, 170, 170)"); // 5, 9, 8, 4, 5,    грань 6
	ind++;
	colors[ind] = new THREE.Color("rgb(200, 200, 200)"); // 6, 10, 9, 5, 6,   грань 7
	ind++;
	colors[ind] = new THREE.Color("rgb(170, 170, 170)"); // 7, 11, 10, 6, 7,  грань 8
	ind++;
	
	// Pavilion
	colors[ind] = new THREE.Color("rgb(250, 190, 190)"); // 12, 11, 8, 12,    грань 9
	ind++;
	colors[ind] = new THREE.Color("rgb(190, 250, 190)"); // 12, 8, 9, 12,     грань 10
	ind++;
	colors[ind] = new THREE.Color("rgb(190, 190, 250)"); // 12, 9, 10, 12,    грань 11
	ind++;
	colors[ind] = new THREE.Color("rgb(150, 250, 250)"); // 12, 10, 11, 12,   грань 12
};