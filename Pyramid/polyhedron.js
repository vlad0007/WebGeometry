// Каждая грань (фасета) трехмерной
// модели представлена нижеприведенной структурой Polygon.
// Структура представляет собой часть плоскости,
// ограниченной многоугольником (полигоном).
// Каждый многоугольник определяется набором 
// индексов вершин принадлежащих модели.
// Например, если грань ограничена многоугольником, 
// который можно обойти путем прохода через вершины 8, 16, 23 и 15,
// то эта грань должна быть записана в файле index.js
// в виде 8, 16, 23, 15, 8. Повторение вершины 8 говорит
// о том, что обход данной грани закончен.
// Все грани (полигоны) 3D модели обходим против часовой стрелки.
// По умолчанию в OpenGL, который является прототипом для WebGL, та сторона полигона, 
// для которой порядок перечисления вершин совпадает с обходом против часовой стрелки, 
// считается лицевой. На модель при этом, конечно, смотрим снаружи. 
// Чтобы все стало более понятным следует запустить программы Pyramid.html и Pyramid_text.html
// в которых приведена подробная нумерация вершин модели "пирамида".

function Polygon()
{
	this.IndexFacet = [];	  // индексы вершин грани без дублированнойй первой вершиной грани
	this.IndexFacet_1 = [];   // индексы вершин грани с дублированной первой вершиной грани
	this.VertexFacet;         // массив, содержащий координаты каждой вершины данной грани
	this.EdgeFacet = [];      // массив, содержащий индексы, определяющие ребро грани (каждый элемент массива содержит два индекса вершин)
	this.IndexTriangle = [];  // массив, содержащий индексы каждого треугольника, из которых состоит грань
	this.VertexTriangle = []; // массив, содержащий координаты каждой вершины каждого треугольника, из которых состоит грань
	this.Faces = [];          // смотри комментарий ниже в тексте программы
}
	
function IndexTriangle (ind0, ind1, ind2)
{
	this.ind0 = ind0;
	this.ind1 = ind1;
	this.ind2 = ind2;
}

function VertexTriangle()
{
	this.v1 = [3];
	this.v2 = [3];
	this.v3 = [3];
}

function VertexFacet()
{
	this.vertexes = []; 
	this.color;
}

function Edge()
{
	this.v_begin;
	this.v_end;
};

var plgs = [];  // Массив граней из которых состоит модель.

var faces = []; // Каждый элемент этого массива состоит из трех индексов
                // треугольников, полученных после триангуляции текущей грани. 
                // Количество элементов массива равно общему количеству
                // треугольников полученных в результате триангуляции    
                // всех граней модели. 
				//  !!! В three.js используется объект Geometry и связанный с ним массив faces.
				// В дальнейшем мы используем массив faces следующим образом, например так:
				//    geometry_freshnel.faces = faces;
				//        или так:
				//    geometry_dispersion.faces = faces;  

// Создаем массив граней модели - заполняем массивы plgs и faces.
function CreatePolyhedron()
{
	var vertex = [];
	
	var i, j;
	var el = 0;
	
	// Координаты всех вершин модели по X, по Y и по Z записаны последовательно.
	// Поэтому необходимо брать последовательно по три числа и сформировать из
	// каждой тройки чисел координату вершины модели в виде Point3D.
	for (i = 0; i < vertices.length/3; i++)
	{
		var pt = new Point3D();
		for (j = 0; j < 3; j++)
		{
			pt[0] = vertices[el];
			pt[1] = vertices[el + 1];
			pt[2] = vertices[el + 2];
		}
		vertex.push(pt);
		el = el + 3;
	}
	
	var index;
	var index_begin;
	var i_index = 0; // Номер индекса, проходит по всем вершинам огранки

	var iface = 0;  
	var iPolyg = 0;
	var indElement = 0;
	i = 0;
	
	for (;;) // Цикл по всем полигонам
	{
		// Полигон 
		var plg_out = new Polygon();  

		index = index_cut[i_index];
		i_index++; // Сразу делаем инкремент - вдруг нарвемся на "break" !!!
		if (index == -100)
			break;	// Прошли по всем полигонам

		index_begin = index;
		plg_out.IndexFacet.push(index);
		plg_out.IndexFacet_1.push(index);
		for (;;)
		{
			// В текущем полигоне заполняем массив индексов его вершин
			index = index_cut[i_index];
			plg_out.IndexFacet.push(index);
			if (index != index_begin)
			{
				plg_out.IndexFacet_1.push(index);
			}

			i_index++; // Берем следующую вершину текущей грани

			//  Для каждой вершины, принадлежащей текущему полигону
			// записываем номер полигона, которому прнадлежит вершина полигона.
			// Нумерация начинается с 0.
			//  Например, предположим модель имеет 34 грани. Если нулевая грань имеет 4 вешины,
			// первая грань - 3 вершины, вторая - 5 пять вершин и т.д.
			// то для огранки записываем последовательность 
			//                    0, 0, 0, 0,   // 
			//                    1, 1, 1,
			//                    2, 2, 2, 2, 2,
			//                    ...........
			//                    ...........
			//                    32, 32, 32, 32, 32,
			//                    33, 33, 33,
			
			plg_out.Faces.push(iface);

			if (index == index_begin)
			{
				// Нашли признак конца вершин для текущей грани
				var vrtF = new VertexFacet();
				var k = 0;
				for (k = 0; k < plg_out.IndexFacet.length; k++)
				{
					var x = vertex[plg_out.IndexFacet[k]][0];
					var y = vertex[plg_out.IndexFacet[k]][1];
					var z = vertex[plg_out.IndexFacet[k]][2];
					var pt = new Point3D(x, y, z);
					vrtF.vertexes.push(pt);
					
				}
				plg_out.VertexFacet = vrtF;
				iface = iface + 1;	// подготовка для следующей грани
				break; // все вершины текущей грани прошли
			}
		}
		
		var i_triangle = 0; // текущий треугольник в данном полигоне
		// nTriangles - кол-во треугольников в текущем полигоне
		var nTriangles = plg_out.IndexFacet_1.length - 2; 
		for (i_triangle = 0; i_triangle < nTriangles; i_triangle++)
		{
			// Триангуляция выпуклого многоугольника.
			var ind_begin = plg_out.IndexFacet_1[0];
			var ind_1 = plg_out.IndexFacet_1[i_triangle+1];
			var ind_2 = plg_out.IndexFacet_1[i_triangle+2];
			var index = new IndexTriangle();
			index.ind0 = ind_begin; // индекс первой вершины
			index.ind1 = ind_1      // индекс второй вершины
			index.ind2 = ind_2;     // индекс третьей вершины
			// Здесь единственное место в данном файле используется объект из библиотеки three.js
			// Если вы используете какую-либо другую библиотеку трехмерной графики,
			// то соответствующим образом измените следующие две строчки программы.
			plg_out.IndexTriangle.push(new THREE.Face3(ind_begin, ind_1, ind_2));  // для three.js
			faces.push(new THREE.Face3(ind_begin, ind_1, ind_2));                  // для three.js
			
			// VertexTriangle vrt;
			var vrt = new VertexTriangle();
			vrt.v1[0] = vertex[ind_begin][0];
			vrt.v1[1] = vertex[ind_begin][1];
			vrt.v1[2] = vertex[ind_begin][2];

			vrt.v2[0] = vertex[ind_1][0];
			vrt.v2[1] = vertex[ind_1][1];
			vrt.v2[2] = vertex[ind_1][2];

			vrt.v3[0] = vertex[ind_2][0];
			vrt.v3[1] = vertex[ind_2][1];
			vrt.v3[2] = vertex[ind_2][2];

			plg_out.VertexTriangle.push(vrt);
		}
		
		// Определяем пары индексов огранки, которые задают ребра огранки
		var nIndLines = 0;
		// Цикл по индексам грани 
		for (j = 0; j < plg_out.IndexFacet.length - 1; j++)
		{
			var edge = new Edge();
			edge.v_begin = plg_out.IndexFacet[j];
			if (j < plg_out.IndexFacet.length - 1)
				edge.v_end = plg_out.IndexFacet[j+1];
			else
				edge.v_end = plg_out.IndexFacet[0];
			
			plg_out.EdgeFacet.push(edge);
		}
		indElement = indElement + nTriangles + 2;
		plgs.push(plg_out);
		iPolyg++;
	}  // здесь действие переменной i_index, ранее объявленной, закончилось
}