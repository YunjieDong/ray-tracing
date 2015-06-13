const SPHERE      = 0;
const BOX         = 1;
const GNDPLANE    = 2;
const CYLINDER    = 3;
const TRIANGLE    = 4;
const BLOBBIES    = 5;
//==============================================================================
// Vertex shader program
var VSHADER_SOURCE =
  'attribute vec4 a_Position;\n' +
  'attribute vec4 a_Color;\n' +
  'attribute vec4 a_Normal; \n' +
  'attribute vec2 a_TexCoord;\n' +

  'uniform vec3 u_Kd; \n' +
  'uniform mat4 u_ViewMatrix;\n' +
  'uniform mat4 u_ProjMatrix;\n' +
  'uniform mat4 u_ModelMatrix;\n' +
  'uniform mat4 u_NormalMatrix; \n' +
  'uniform int u_flag;\n' +

  'varying vec2 v_TexCoord;\n' +
  'varying vec4 v_Position; \n' +
  'varying vec3 v_Normal; \n' +
  'varying vec4 v_Color;\n' +
  'varying vec3 v_Kd; \n' +
  'void main() {\n' +
  '  if(u_flag > 0) {  \n'+
  '      gl_Position = a_Position;\n' +
  '  } \n'+
  '  else { \n'+
  '      gl_Position = u_ProjMatrix * u_ViewMatrix * u_ModelMatrix * a_Position; \n'+
  '      v_Position = u_ModelMatrix * a_Position; \n' +
  '      v_Normal = normalize(vec3(u_NormalMatrix * a_Normal));\n' +
  '  } \n' +
  '	 v_Kd = u_Kd; \n' +
  '  v_TexCoord = a_TexCoord;\n' +
  '  v_Color = a_Color;\n' +
  '}\n';

// Fragment shader program
var FSHADER_SOURCE =
  '#ifdef GL_ES\n' +
  'precision mediump float;\n' +
  '#endif\n' +

  // first light source------------------------------------------
  'uniform vec4 u_Lamp0Pos;\n' + 			// Phong Illum: position
  'uniform vec3 u_Lamp0Amb;\n' +   		    // Phong Illum: ambient
  'uniform vec3 u_Lamp0Diff;\n' +           // Phong Illum: diffuse
  'uniform vec3 u_Lamp0Spec;\n' +			// Phong Illum: specular
  // second light source------------------------------------------
  'uniform vec4 u_Lamp1Pos;\n' + 			// Phong Illum: position
  'uniform vec3 u_Lamp1Amb;\n' +   		    // Phong Illum: ambient
  'uniform vec3 u_Lamp1Diff;\n' +           // Phong Illum: diffuse
  'uniform vec3 u_Lamp1Spec;\n' +			// Phong Illum: specular

  // first material definition
  'uniform vec3 u_Ke;\n' +					// Phong Reflectance: emissive
  'uniform vec3 u_Ka;\n' +					// Phong Reflectance: ambient
  'uniform vec3 u_Ks;\n' +				    // Phong Reflectance: specular

  'uniform vec4 u_eyePosWorld; \n' + 		// Camera/eye location in world coords.
  'uniform int u_isTexture; \n' +			// texture/not-texture flag
  'uniform int u_isLighting; \n' +
  'uniform int u_Light0; \n' +
  'uniform int u_Light1; \n' +
  'uniform sampler2D u_Sampler;\n' +	    // our 2D texture-addr-maker

  'varying vec2 v_TexCoord;\n' +
  'varying vec4 v_Color;\n' +
  'varying vec3 v_Normal;\n' +				// Find 3D surface normal at each pix
  'varying vec4 v_Position;\n' +			// pixel's 3D pos too -- in 'world' coords
  'varying vec3 v_Kd;	\n' +				// Find diffuse reflectance K_d per pix

  'void main() {\n' +
  '  vec3 normal = normalize(v_Normal); \n' +
  '  vec3 lightDirection = normalize(u_Lamp0Pos.xyz - v_Position.xyz);\n' +
  '  vec3 lightDirection1 = normalize(u_Lamp1Pos.xyz - v_Position.xyz);\n' +
  '  float nDotL = max(dot(lightDirection, normal), 0.0); \n' +
  '  float nDotL1 = max(dot(lightDirection1, normal), 0.0); \n' +
  '  vec3 eyeDirection = normalize(u_eyePosWorld.xyz - v_Position.xyz); \n' +
  '  vec3 H = normalize(lightDirection + eyeDirection); \n' +
  '  vec3 H1 = normalize(lightDirection1 + eyeDirection); \n' +
  '  float nDotH = max(dot(H, normal), 0.0); \n' +
  '  float nDotH1 = max(dot(H1, normal), 0.0); \n' +

  '  float e02 = nDotH*nDotH; \n' +
  '  float e04 = e02*e02; \n' +
  '  float e08 = e04*e04; \n' +
  '	 float e16 = e08*e08; \n' +
  '	 float e32 = e16*e16; \n' +
  '	 float e64 = e32*e32; \n' +
  //------------------------------------------
  '  float e021 = nDotH1*nDotH1; \n' +
  '  float e041 = e021*e021; \n' +
  '  float e081 = e041*e041; \n' +
  '	 float e161 = e081*e081; \n' +
  '	 float e321 = e161*e161; \n' +
  '	 float e641 = e321*e321; \n' +
  // Calculate the final color from diffuse reflection and ambient reflection
  '	 vec3 emissive = u_Ke;' +
  '  vec3 ambient = u_Lamp0Amb * u_Ka;\n' +
  '  vec3 diffuse = u_Lamp0Diff * v_Kd * nDotL;\n' +
  '	 vec3 speculr = u_Lamp0Spec * u_Ks * e64 * e64;\n' +
  '  vec3 sum = vec3(emissive + ambient + diffuse + speculr);\n' +

  //------------------------------------------
  '	 vec3 emissive1 = u_Ke;' +
  '  vec3 ambient1 = u_Lamp1Amb * u_Ka;\n' +
  '  vec3 diffuse1 = u_Lamp1Diff * v_Kd * nDotL1;\n' +
  '	 vec3 speculr1 = u_Lamp1Spec * u_Ks * e641 * e641;\n' +
  '  vec3 sum1 = emissive1 + ambient1 + diffuse1 + speculr1;\n' +

  '  if(u_isTexture > 0) {  \n' +				// pixel color comes from texture-map,
  '  	 gl_FragColor = texture2D(u_Sampler, v_TexCoord);\n' +
  '  } \n' +
  '  else { \n' +
  '      if(u_isLighting > 0) {  \n' +
  '          if(u_Light0 > 0) {  \n' +
  '              if(u_Light1 > 0) {  \n' +
  '	 	             gl_FragColor = vec4(sum + sum1, 1.0); \n' +
  '              } \n' +
  '              else { \n' +
  '                  gl_FragColor = vec4(sum, 1.0); \n' +
  '              } \n' +
  '          } \n' +
  '          else { \n' +
  '              if(u_Light1 > 0) {  \n' +
  '                  gl_FragColor = vec4(sum1, 1.0); \n'+
  '              } \n' +
  '              else { \n' +
  '                  gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0); \n' +
  '              } \n' +
  '          } \n' +
  '      } \n' +
  '      else { \n' +
  '          gl_FragColor = v_Color; \n' +
  '      } \n' +
  '  } \n' +
  '}\n';

//Global Vars:---------------------------------------------------------------
//-----------'Uniform' values (sent to the GPU)
var u_flag = 0;
var u_isLighting = 1;
var light = 1;
var u_isTexture = 0;					// ==0 false--use fixed colors in frag shader
                                        // ==1 true --use texture-mapping in frag shader
var u_isTextureID = 0;			  // GPU location of this uniform var
var u_Light0 = 1;
var u_Light1 = 1;
var light0 = 1;
var light1 = 1;
var shadowFlag = 1;
var scene = 1;

//-----------Global vars for mouse click-and-drag for rotation.
var isDrag=false;		// mouse-drag: true when user holds down mouse button
var xMclik=0.0;			// last mouse button-down position (in CVV coords)
var yMclik=0.0;
var xMdragTot=0.0;	// total (accumulated) mouse-drag amounts (in CVV coords).
var yMdragTot=0.0;

var qNew = new Quaternion(0,0,0,1); // most-recent mouse drag's rotation
var qTot = new Quaternion(0,0,0,1);	// 'current' orientation (made from qNew)
var quatMatrix = new Matrix4();				// rotation matrix, made from latest qTot

//-----------Ray Tracer Objects:
var myPic = new CImgBuf(256,256);	// create RGB image buffer object, and
var isTracing = 1;									// ==1 when tracing; else 0

//------------------------
var floatsPerVertex = 11;
var xcount = 101;			// # of lines to draw in x,y to make the grid.
var ycount = 101;
var xymax = 50.0;			// grid size; extends to cover +/-xymax in x and y.

var g_EyeX = 0.0, g_EyeY = 2.5, g_EyeZ = 1.0;
var g_AtX = 0.0, g_AtY = 1.0, g_AtZ = -5.0;
var theta = 0;
var zfar = 100;

var lampX = 0.0;
var lampY = 26.0;
var lampZ = -8.5;

//==============================================================================
CRay = function(){
    this.orig = vec4.fromValues(g_EyeX, g_EyeY, g_EyeZ, 1.0);
    this.dir = vec4.fromValues(0,0,-1.0,0);
}
CRay.prototype.constructor = CRay;

myRay = new CRay();
//==============================================================================
CHit = function(){
    this.hitItem = -1;
    this.th = zfar;
    this.modelHitPt = vec4.create();        // the 'hit point' expressed in model coords.
    this.colr = vec3.create();
    this.isEntering = true;     // true iff ray origin was OUTSIDE the hitItem.
    
    this.hitPt = vec4.fromValues(zfar, zfar, zfar, 1.0);    // World-space location where the ray pierced
    this.surfNorm = vec4.create();      // World-space surface-normal vector at the hit point
    this.viewN = vec4.create();     // Unit-length vector from hitPt back towards the origin of the ray
    this.depth = 0;
}
CHit.prototype.constructor = CHit();
//==============================================================================
CHitList = function(){
    this.pierce = new Array();
    this.iEnd = 1;
    this.iNearest = 0;
    this.isShadowRay = false;
}
CHitList.prototype.constructor = CHitList;
CHitList.prototype.initList = function(bkColr){
    this.pierce[0] = new CHit();
    vec4.copy(this.pierce[0].colr, bkColr);
}
CHitList.prototype.initShadow = function(){
    this.isShadowRay = true;
    this.pierce[0].th = -1;
}
CHitList.prototype.findNearestHit = function(){
    switch(this.iEnd){
        case 0:
            break;
        case 1:
            this.iNearest = 0;
            break;
        case 2:
            this.iNearest = 1;
            break;
        case 3:
            if (this.pierce[1].th < this.pierce[2].th){
                this.iNearest = 1;
            }
            else{
                this.iNearest = 2;
            }
            break;
        default:
            this.iNearest = 1;
            for (var i=2; i<this.iEnd; i++){
                if (this.pierce[i].th < this.pierce[this.iNearest].th){
                    this.iNearest = i;
                }
            }
            break;
    }
}

myHitList = new CHitList();
//==============================================================================
CCamera = function(fov, aspect, near, far){
    this.projectMatrix = mat4.create();
    mat4.perspective(this.projectMatrix, fov, aspect, near, far);
    var left = vec4.fromValues(-1, 0, -1, 1);
    var right = vec4.fromValues(1, 0, -1, 1);
    var top = vec4.fromValues(0, 1, -1, 1);
    var bottom = vec4.fromValues(0, -1, -1, 1);
    
    // create and calculate the inverse matrix from the project matrix
    var invProjMatrix = mat4.create();
    mat4.invert(invProjMatrix, this.projectMatrix);
    vec4.transformMat4(left, left, invProjMatrix);
    vec4.transformMat4(right, right, invProjMatrix);
    vec4.transformMat4(top, top, invProjMatrix);
    vec4.transformMat4(bottom, bottom, invProjMatrix);
    // keep w = 1
    vec4.scale(left, left, 1.0 / left[3]);
    vec4.scale(right, right, 1.0 / right[3]);
    vec4.scale(top, top, 1.0 / top[3]);
    vec4.scale(bottom, bottom, 1.0 / bottom[3]);
    
    this.iLeft = left[0];
    this.iRight = right[0];
    this.iBot = bottom[1];
    this.iTop = top[1];
    this.znear = -1.0;
    
    this.xmax = 256;
    this.ymax = 256;
    this.ufrac = (this.iRight-this.iLeft)/this.xmax;
    this.vfrac = (this.iTop-this.iBot)/this.ymax;
    
    this.vrp = vec4.create();
    this.lookAtPt = vec4.create();
    this.vup = vec4.fromValues(0.0, 1.0, 0.0, 0.0);
    
    this.N = vec4.create();
    this.U = vec4.create();
    this.V = vec4.create();
    this.LL = vec4.create();
    
//    this.LL = vec4.fromValues(this.iLeft, this.iBot, this.znear, 1.0);
    this.p0 = vec4.fromValues(this.iLeft+0.5*this.ufrac, this.iBot+0.5*this.vfrac, this.znear, 1.0);
    this.uvnPix = vec4.fromValues(0,0,0,1.0);

}
CCamera.prototype.constructor = CCamera;
CCamera.prototype.setFrustum = function(){
    vec4.copy(this.vrp, myRay.orig);
    this.lookAtPt = vec4.fromValues(g_AtX, g_AtY, g_AtZ, 1.0);
    
    this.N = vec4.sub(vec4.create(), this.vrp, this.lookAtPt);
    vec4.normalize(this.N,this.N);
    
    this.U = vec3.cross(vec3.create(), this.vup, this.N);
    this.U = vec4.fromValues(this.U[0], this.U[1], this.U[2], 0.0);
    vec4.normalize(this.U,this.U);
    
    this.V = vec3.cross(vec3.create(), this.N, this.U);
    this.V = vec4.fromValues(this.V[0], this.V[1], this.V[2], 0.0);
    vec4.normalize(this.V,this.V);

    vec4.copy(this.LL, this.vrp);
    vec4.scaleAndAdd(this.LL, this.LL, this.U, this.iLeft);
    vec4.scaleAndAdd(this.LL, this.LL, this.V, this.iBot);
    vec4.scaleAndAdd(this.LL, this.LL, this.N, this.znear);
}
CCamera.prototype.makeEyeRay = function(eyeRay){
    var xyzPos = vec4.create();
    
    vec4.scaleAndAdd(xyzPos, xyzPos, this.U, this.uvnPix[0]);
    vec4.scaleAndAdd(xyzPos, xyzPos, this.V, this.uvnPix[1]);
    vec4.scaleAndAdd(xyzPos, xyzPos, this.N, this.znear);
    
    vec4.copy(eyeRay.dir, xyzPos);
}

//==============================================================================
CGeom = function(type){
    this.enabled = true;
    this.shapeType = type;
    this.matlNum = -1;
    this.shininess = 0.0;
    this.transparency = 0.0;
    this.refractIndex = 1.0;
    
    //Grid------------------------------------------
    this.zVal = -5.0;
    this.xgap = 2*xymax/(xcount-1);
    this.ygap = 2*xymax/(ycount-1);
    this.linewidth = 0.1;
    this.lineColr = vec3.fromValues(0.0, 0.0, 0.0);
    this.gapColr = vec3.fromValues(1.0, 1.0, 1.0);
    
    //Sphere------------------------------------------
    this.center = vec4.fromValues(0,0,0,1);
    this.radius = 1.0;
    this.sqrRadius = this.radius * this.radius;
    this.sphereColr = vec3.fromValues(1.0, 1.0, 0.5);
    
    //Box------------------------------------------
    this.boxCenter = vec4.fromValues(0,0,0,1);
    this.length = 1.0;
    this.hLength = this.length/2;
    this.boxColr = vec3.fromValues(0.94, 0.27, 0.31);
    this.dist = new Array();
    this.dist[0] = this.boxCenter[0] - this.hLength;
    this.dist[1] = this.boxCenter[0] + this.hLength;
    this.dist[2] = this.boxCenter[1] - this.hLength;
    this.dist[3] = this.boxCenter[1] + this.hLength;
    this.dist[4] = this.boxCenter[2] - this.hLength;
    this.dist[5] = this.boxCenter[2] + this.hLength;
    
    //Cylinder------------------------------------------
    this.Zmin = 0;
    this.Zmax = 1;
    this.topRadius = 1.0;
    this.cylColr = vec3.fromValues(249/255, 205/255, 173/255);
    
    this.th = 0;
    this.pos = vec4.fromValues(0,0,0,1.0);
    this.colr = vec3.fromValues(0,0,0);
    
    this.model2world = mat4.create();
    this.world2model = mat4.create();
    mat4.invert(this.world2model, this.model2world);
    this.w2mTranspose = mat4.create();
}
CGeom.prototype.constructor = CGeom;
CGeom.prototype.findNormal = function(HitObj){
    mat4.transpose(this.w2mTranspose, this.world2model);
    switch(this.shapeType){
        case GNDPLANE:
            HitObj.surfNorm = vec4.fromValues(0,0,1,0);
            vec4.transformMat4(HitObj.surfNorm, HitObj.surfNorm, this.w2mTranspose);
            HitObj.surfNorm[3] = 0.0;
            vec4.normalize(HitObj.surfNorm, HitObj.surfNorm);
            break;
        case SPHERE:
            vec4.sub(HitObj.surfNorm, HitObj.modelHitPt, this.center);
            vec4.transformMat4(HitObj.surfNorm, HitObj.surfNorm, this.w2mTranspose);
            HitObj.surfNorm[3] = 0.0;
            vec4.normalize(HitObj.surfNorm, HitObj.surfNorm);
            break;
        case BOX:
            if (HitObj.modelHitPt[0] == this.dist[0]){
                HitObj.surfNorm = vec4.fromValues(-1,0,0,0);
            }
            else if (HitObj.modelHitPt[0] == this.dist[1]){
                HitObj.surfNorm = vec4.fromValues(1,0,0,0);
            }
            else if (HitObj.modelHitPt[1] == this.dist[2]){
                HitObj.surfNorm = vec4.fromValues(0,-1,0,0);
            }
            else if (HitObj.modelHitPt[1] == this.dist[3]){
                HitObj.surfNorm = vec4.fromValues(0,1,0,0);
            }
            else if (HitObj.modelHitPt[2] == this.dist[4]){
                HitObj.surfNorm = vec4.fromValues(0,0,-1,0);
            }
            else if (HitObj.modelHitPt[2] == this.dist[5]){
                HitObj.surfNorm = vec4.fromValues(0,0,1,0);
            }
            vec4.transformMat4(HitObj.surfNorm, HitObj.surfNorm, this.w2mTranspose);
            HitObj.surfNorm[3] = 0.0;
            vec4.normalize(HitObj.surfNorm, HitObj.surfNorm);
            break;
        case CYLINDER:
            if (HitObj.modelHitPt[2] == this.Zmin){
                HitObj.surfNorm = vec4.fromValues(0,0,-1,0);
            }
            else if (HitObj.modelHitPt[2] == this.Zmax){
                HitObj.surfNorm = vec4.fromValues(0,0,1,0);
            }
            else{
                var sm = this.topRadius - 1;
                HitObj.surfNorm = vec4.fromValues(HitObj.modelHitPt[0], HitObj.modelHitPt[1], -sm * (1+sm*HitObj.modelHitPt[2]), 0);
            }
            vec4.transformMat4(HitObj.surfNorm, HitObj.surfNorm, this.w2mTranspose);
            HitObj.surfNorm[3] = 0.0;
            vec4.normalize(HitObj.surfNorm, HitObj.surfNorm);
            break;
            
    }
}
CGeom.prototype.drawGrid = function(){
    makeGroundGrid();
}
CGeom.prototype.traceGrid = function(myRay, myList, index){
// return -1 if ray MISSES the plane
// return  0 if ray hits BETWEEN lines
// return  1 if ray hits ON the lines
    var myHitGrid = new CHit();
    var torg = vec4.transformMat4(vec4.create(), myRay.orig, this.world2model);
    var tdir = vec4.transformMat4(vec4.create(), myRay.dir, this.world2model);
    vec4.scale(torg, torg, 1.0 / torg[3]);
    
    this.th = (0 - torg[2])/tdir[2];
    if (this.th < 0.0){
        return -1;      //miss the plane
    }
    this.pos[0] = torg[0] + tdir[0]*this.th;
    this.pos[1] = torg[1] + tdir[1]*this.th;
    this.pos[2] = torg[2] + tdir[2]*this.th;
    
    myList.pierce[myList.iEnd] = myHitGrid;
    myList.iEnd += 1;
    myHitGrid.th = this.th;
    myHitGrid.hitItem = index;
    myHitGrid.isEntering = true;
    if (myList.isShadowRay == false){
        vec4.scale(myHitGrid.viewN, myRay.dir, -1.0);       // (store reversed ray direction)
        vec4.normalize(myHitGrid.viewN, myHitGrid.viewN);
        vec4.copy(myHitGrid.modelHitPt, this.pos);
    }
    
    var xfrac = this.pos[0] - Math.floor(this.pos[0]/this.xgap);
    var yfrac = this.pos[1] - Math.floor(this.pos[1]/this.ygap);
    if (xfrac < this.linewidth || yfrac < this.linewidth){
        myHitGrid.colr = this.lineColr;
    }
    else{
        myHitGrid.colr = this.gapColr;
    }
}
CGeom.prototype.drawSphere = function(sColr){
    makeSphere(sColr);
}
CGeom.prototype.traceSphere = function(myRay, myList, index){
// return -1 if ray MISSES the sphere
// return  0 if ray hits
    var myHitSphere = new CHit();
    var myHitSphere1 = new CHit();
    var torg = vec4.transformMat4(vec4.create(), myRay.orig, this.world2model);
    var tdir = vec4.transformMat4(vec4.create(), myRay.dir, this.world2model);
    vec4.scale(torg, torg, 1.0 / torg[3]);
    
    var DL2 = vec4.squaredLength(tdir);     //find tdir squared length, DL2
    var r2s = vec4.create();        // r2s = vector from ray origin to sphere center
    vec4.sub(r2s, this.center, torg);
    var L2 = vec4.squaredLength(r2s);       // Find squared length of that vector, L2
    var tcaS = vec4.dot(tdir, r2s);     // SCALED tca; tca*length(tdir).
    
    if (L2 > this.sqrRadius){       // the ray begins OUTSIDE the sphere ( because L2 > radius^2)
        if (tcaS < 0.0){        // MISS! sphere is BEHIND the ray! No hit points.
            return -1;
        }
        var tca2 = tcaS*tcaS / DL2;
        var LM2 = L2 - tca2;        //find LM2; the squared distance from sphere center to chord mid-point
        if (LM2 > this.sqrRadius){      // MISS! ray point nearest to sphere center > sphere radius!
            return -1;
        }
        if (myList.isShadowRay == true){     //we don't care about WHERE we hit
            myList.pierce[myList.iEnd] = myHitSphere;
            myList.iEnd += 1;
            myHitSphere.th = 1.0;     // set bogus hit distance
            myHitSphere.hitItem = index;
            return;
        }
        var Lhc2 = this.sqrRadius - LM2;        // SQUARED half-chord length.
        //  Record first of two hit points;-----------------------------------------
        this.th = (tcaS/DL2) - Math.sqrt(Lhc2/DL2);     //at smaller t we ENTER the sphere
        this.pos[0] = torg[0] + tdir[0]*this.th;
        this.pos[1] = torg[1] + tdir[1]*this.th;
        this.pos[2] = torg[2] + tdir[2]*this.th;
        
        myList.pierce[myList.iEnd] = myHitSphere;
        myList.iEnd += 1;
        myHitSphere.th = this.th;
        myHitSphere.colr = this.sphereColr;
        myHitSphere.isEntering = true;
        myHitSphere.hitItem = index;
        vec4.scale(myHitSphere.viewN, myRay.dir, -1.0);
        vec4.normalize(myHitSphere.viewN, myHitSphere.viewN);
        vec4.copy(myHitSphere.modelHitPt, this.pos);
        
        //  Record second of two hit points:----------------------------------------
        this.th = (tcaS/DL2) + Math.sqrt(Lhc2/DL2);     //at larger t we LEAVE the sphere
        this.pos[0] = torg[0] + tdir[0]*this.th;
        this.pos[1] = torg[1] + tdir[1]*this.th;
        this.pos[2] = torg[2] + tdir[2]*this.th;
        
        myList.pierce[myList.iEnd] = myHitSphere1;
        myList.iEnd += 1;
        myHitSphere1.th = this.th;
        myHitSphere1.colr = this.sphereColr;
        myHitSphere1.isEntering = false;
        myHitSphere1.hitItem = index;
        vec4.scale(myHitSphere1.viewN, myRay.dir, -1.0);
        vec4.normalize(myHitSphere1.viewN, myHitSphere1.viewN);
        vec4.copy(myHitSphere1.modelHitPt, this.pos);
    }
    else{
        var tca2 = tcaS*tcaS / DL2;
        var LM2 = L2 - tca2;
        if (LM2 > this.sqrRadius){      // MISS! ray point nearest to sphere center > sphere radius!
            return -1;
        }
        var Lhc2 = this.sqrRadius - LM2;
        this.th = (tcaS/DL2) + Math.sqrt(Lhc2/DL2);
        this.pos[0] = torg[0] + tdir[0]*this.th;
        this.pos[1] = torg[1] + tdir[1]*this.th;
        this.pos[2] = torg[2] + tdir[2]*this.th;
        
        myList.pierce[myList.iEnd] = myHitSphere;
        myList.iEnd += 1;
        myHitSphere.th = this.th;
        myHitSphere.colr = this.sphereColr;
        myHitSphere.isEntering = false;
        myHitSphere.hitItem = index;
        vec4.scale(myHitSphere.viewN, myRay.dir, -1.0);
        vec4.normalize(myHitSphere.viewN, myHitSphere.viewN);
        vec4.copy(myHitSphere.modelHitPt, this.pos);
    }
    
}
CGeom.prototype.drawBox = function(){
    makeCube();
}
CGeom.prototype.traceBox = function(myRay, myList, index){
    var boxList = new Array();
    var torg = vec4.transformMat4(vec4.create(), myRay.orig, this.world2model);
    var tdir = vec4.transformMat4(vec4.create(), myRay.dir, this.world2model);
    vec4.scale(torg, torg, 1.0 / torg[3]);
    
    for (var i=0; i<2; i++){
        this.th = (this.dist[i] - torg[0])/tdir[0];
        vec4.scale(this.pos, tdir, this.th);
        vec4.add(this.pos, this.pos, torg);
        if (this.pos[1]>=this.dist[2] && this.pos[1]<=this.dist[3] && this.pos[2]>=this.dist[4] && this.pos[2]<=this.dist[5]){
            myList.pierce[myList.iEnd] = new CHit();
            myList.pierce[myList.iEnd].th = this.th;
            myList.pierce[myList.iEnd].hitItem = index;
            myList.pierce[myList.iEnd].colr = this.boxColr;
            vec4.scale(myList.pierce[myList.iEnd].viewN, myRay.dir, -1.0);     // (store reversed ray direction)
            vec4.normalize(myList.pierce[myList.iEnd].viewN, myList.pierce[myList.iEnd].viewN);
            vec4.copy(myList.pierce[myList.iEnd].modelHitPt, this.pos);
            boxList.push(myList.pierce[myList.iEnd]);
            myList.iEnd += 1;
        }
    }
    for (var i=0; i<2; i++){
        this.th = (this.dist[i+2] - torg[1])/tdir[1];
        vec4.scale(this.pos, tdir, this.th);
        vec4.add(this.pos, this.pos, torg);
        if (this.pos[0]>=this.dist[0] && this.pos[0]<=this.dist[1] && this.pos[2]>=this.dist[4] && this.pos[2]<=this.dist[5]){
            myList.pierce[myList.iEnd] = new CHit();
            myList.pierce[myList.iEnd].th = this.th;
            myList.pierce[myList.iEnd].hitItem = index;
            myList.pierce[myList.iEnd].colr = this.boxColr;
            vec4.scale(myList.pierce[myList.iEnd].viewN, myRay.dir, -1.0);     // (store reversed ray direction)
            vec4.normalize(myList.pierce[myList.iEnd].viewN, myList.pierce[myList.iEnd].viewN);
            vec4.copy(myList.pierce[myList.iEnd].modelHitPt, this.pos);
            boxList.push(myList.pierce[myList.iEnd]);
            myList.iEnd += 1;
        }
    }
    for (var i=0; i<2; i++){
        this.th = (this.dist[i+4] - torg[2])/tdir[2];
        vec4.scale(this.pos, tdir, this.th);
        vec4.add(this.pos, this.pos, torg);
        if (this.pos[0]>=this.dist[0] && this.pos[0]<=this.dist[1] && this.pos[1]>=this.dist[2] && this.pos[1]<=this.dist[3]){
            myList.pierce[myList.iEnd] = new CHit();
            myList.pierce[myList.iEnd].th = this.th;
            myList.pierce[myList.iEnd].hitItem = index;
            myList.pierce[myList.iEnd].colr = this.boxColr;
            vec4.scale(myList.pierce[myList.iEnd].viewN, myRay.dir, -1.0);     // (store reversed ray direction)
            vec4.normalize(myList.pierce[myList.iEnd].viewN, myList.pierce[myList.iEnd].viewN);
            vec4.copy(myList.pierce[myList.iEnd].modelHitPt, this.pos);
            boxList.push(myList.pierce[myList.iEnd]);
            myList.iEnd += 1;
        }
    }
    if (boxList.length == 1){
        boxList[0].isEntering = false;
    }
    else if (boxList.length == 2){
        if (boxList[0].th > boxList[1].th){
            boxList[0].isEntering = false;
            boxList[1].isEntering = true;
        }
        else{
            boxList[0].isEntering = true;
            boxList[1].isEntering = false;
        }
    }

}
CGeom.prototype.drawCylinder = function(tR){
    makeCylinder(tR);
}
CGeom.prototype.traceCylinder = function(myRay, myList, index){
    var torg = vec4.transformMat4(vec4.create(), myRay.orig, this.world2model);
    var tdir = vec4.transformMat4(vec4.create(), myRay.dir, this.world2model);
    vec4.scale(torg, torg, 1.0 / torg[3]);
    
    var d = (this.topRadius-1) * tdir[2];
    var F = 1 + (this.topRadius-1)*torg[2];
    
    var A = tdir[0]*tdir[0] + tdir[1]*tdir[1] - d*d;
    var B = torg[0]*tdir[0] + torg[1]*tdir[1] - F*d;
    var C = torg[0]*torg[0] + torg[1]*torg[1] - F*F;
    var discriminant = B*B - A*C;
    
    if (discriminant > 0){
        var root = Math.pow(discriminant, 0.5);
        this.th = (-B - root) / A;      //earlier hit
        vec4.scale(this.pos, tdir, this.th);
        vec4.add(this.pos, this.pos, torg);
        if (this.th > 0.0001 && this.pos[2] >= this.Zmin && this.pos[2] <= this.Zmax){
            myList.pierce[myList.iEnd] = new CHit();
            myList.pierce[myList.iEnd].th = this.th;
            myList.pierce[myList.iEnd].hitItem = index;
            myList.pierce[myList.iEnd].colr = this.cylColr;
            myList.pierce[myList.iEnd].isEntering = true;
            vec4.scale(myList.pierce[myList.iEnd].viewN, myRay.dir, -1.0);     // (store reversed ray direction)
            vec4.normalize(myList.pierce[myList.iEnd].viewN, myList.pierce[myList.iEnd].viewN);
            vec4.copy(myList.pierce[myList.iEnd].modelHitPt, this.pos);
            //cylList.push(myList.pierce[myList.iEnd]);
            myList.iEnd += 1;
        }
        this.th = (-B + root) / A;      //second hit
        vec4.scaleAndAdd(this.pos, torg, tdir, this.th);
        if (this.th > 0.0001 && this.pos[2] >= this.Zmin && this.pos[2] <= this.Zmax){
            myList.pierce[myList.iEnd] = new CHit();
            myList.pierce[myList.iEnd].th = this.th;
            myList.pierce[myList.iEnd].hitItem = index;
            myList.pierce[myList.iEnd].colr = this.cylColr;
            myList.pierce[myList.iEnd].isEntering = false;
            vec4.scale(myList.pierce[myList.iEnd].viewN, myRay.dir, -1.0);     // (store reversed ray direction)
            vec4.normalize(myList.pierce[myList.iEnd].viewN, myList.pierce[myList.iEnd].viewN);
            vec4.copy(myList.pierce[myList.iEnd].modelHitPt, this.pos);
            //cylList.push(myList.pierce[myList.iEnd]);
            myList.iEnd += 1;
        }
    }
    
    this.th = (this.Zmin - torg[2]) / tdir[2];      // bottom plane
    vec4.scaleAndAdd(this.pos, torg, tdir, this.th);
    var area = this.pos[0]*this.pos[0] + this.pos[1]*this.pos[1];
    if (this.th > 0.0001 && area < 1){
        myList.pierce[myList.iEnd] = new CHit();
        myList.pierce[myList.iEnd].th = this.th;
        myList.pierce[myList.iEnd].hitItem = index;
        myList.pierce[myList.iEnd].colr = this.cylColr;
        vec4.scale(myList.pierce[myList.iEnd].viewN, myRay.dir, -1.0);     // (store reversed ray direction)
        vec4.normalize(myList.pierce[myList.iEnd].viewN, myList.pierce[myList.iEnd].viewN);
        vec4.copy(myList.pierce[myList.iEnd].modelHitPt, this.pos);
        //cylList.push(myList.pierce[myList.iEnd]);
        myList.iEnd += 1;
    }
    
    this.th = (this.Zmax - torg[2]) / tdir[2];      // top plane
    vec4.scaleAndAdd(this.pos, torg, tdir, this.th);
    area = this.pos[0]*this.pos[0] + this.pos[1]*this.pos[1];
    if (this.th > 0.0001 && area < this.topRadius*this.topRadius){
        myList.pierce[myList.iEnd] = new CHit();
        myList.pierce[myList.iEnd].th = this.th;
        myList.pierce[myList.iEnd].hitItem = index;
        myList.pierce[myList.iEnd].colr = this.cylColr;
        vec4.scale(myList.pierce[myList.iEnd].viewN, myRay.dir, -1.0);     // (store reversed ray direction)
        vec4.normalize(myList.pierce[myList.iEnd].viewN, myList.pierce[myList.iEnd].viewN);
        vec4.copy(myList.pierce[myList.iEnd].modelHitPt, this.pos);
        //cylList.push(myList.pierce[myList.iEnd]);
        myList.iEnd += 1;
    }
    
}

myGrid = new CGeom(GNDPLANE);
mySphere = new CGeom(SPHERE);
mySphereB = new CGeom(SPHERE);
mySphereC = new CGeom(SPHERE);
mySphereD = new CGeom(SPHERE);
mySphereE = new CGeom(SPHERE);
mySphereF = new CGeom(SPHERE);
mySphereG = new CGeom(SPHERE);
mySphereH = new CGeom(SPHERE);
myBox = new CGeom(BOX);
myCylinder = new CGeom(CYLINDER);
myCone = new CGeom(CYLINDER);
myShade = new CGeom(CYLINDER);
myPillar = new CGeom(CYLINDER);
myBottom = new CGeom(BOX);
myCube = new CGeom(BOX);


//==============================================================================
CMatl = function(){
    this.enabled = true;
    this.matlFlags = 0;
    this.kemi = vec3.create();
    this.kamb = vec3.create();
    this.kdif = vec3.create();
    this.kspec = vec3.create();
    this.kmir = vec3.create();
    this.kshiny = 0.0;
    this.kglass = vec3.create();
    this.keta = 0.0;
    this.kdense = 0.0;
    this.solidID = 0;
    this.matlBID = 0;
}
CMatl.prototype.constructor = CMatl;
CMatl.prototype.lerp = function(blend, src0, src1){
    var K_amb = vec4.lerp(vec4.create(), src0.kamb, src1.kamb, blend);
    var K_dif = vec4.lerp(vec4.create(), src0.kdif, src1.kdif, blend);
    var K_spec = vec4.lerp(vec4.create(), src0.kspec, src1.kspec, blend);
    var K_mir = vec4.lerp(vec4.create(), src0.kmir, src1.kmir, blend);
    var K_glass = vec4.lerp(vec4.create(), src0.kglass, src1.kglass, blend);
    
    var temp = (1.0-blend) * src0.keta + blend * src1.keta;
    var K_eta = temp;
    temp = (1.0-blend) * src0.kdense + blend * src1.kdense;
    var K_dense = temp;
    temp = (1.0-blend) * src0.kshiny + blend * src1.kshiny;
    var K_shiny = temp;
}

myMatl = new CMatl();
myMatlB = new CMatl();
myMatlC = new CMatl();
myMatlD = new CMatl();
myMatlE = new CMatl();
myMatlF = new CMatl();
myMatlG = new CMatl();
myMatlH = new CMatl();
myMatlI = new CMatl();
myMatlJ = new CMatl();
myMatlK = new CMatl();
myMatlL = new CMatl();
myMatlM = new CMatl();
myMatlN = new CMatl();
myMatlO = new CMatl();
myMatlP = new CMatl();
//==============================================================================
CLight = function(){
    this.enabled = true;
    this.lightType = 0;
    this.lightPos = vec4.create();
    this.I_a = vec3.create();
    this.I_d = vec3.create();
    this.I_s = vec3.create();
}
CLight.prototype.constructor = CLight;

myLight = new CLight();
HeadLight = new CLight();
//==============================================================================
CScene = function(){
    this.rayCam = myCam;
    this.rayNow = myRay;
    this.rayHits = myHitList;
    
    this.item = new Array();
    this.itemCount = 0;
    
    this.matList = new Array();
    this.matCount = 0;
    
    this.lightList = new Array();
    this.lampCount = 0;
    
    this.bkgndColr = vec3.fromValues(0.2, 0.2, 0.3);
    this.blackColr = vec3.fromValues(0.0, 0.5, 0.0);
    this.depthMax = 0;
    this.isAA = false;
    this.xSuperAA = 4;
    this.ySuperAA = 4;
    
    this.rayHits.initList(this.bkgndColr);
    
}
CScene.prototype.constructor = CScene;
CScene.prototype.addGeom = function(shape){
    this.item.push(shape);
    this.itemCount += 1;
}
CScene.prototype.addMatl = function(matl){
    this.matList.push(matl);
    this.matCount += 1;
}
CScene.prototype.addLamp = function(lamp){
    this.lightList.push(lamp);
    this.lampCount += 1;
}
CScene.prototype.trace = function(rayNow, rayHits, depth){
    for (var i=0; i<this.itemCount; i++){
        switch(this.item[i].shapeType){
            case GNDPLANE:
                this.item[i].traceGrid(rayNow, rayHits, i);
                break;
            case SPHERE:
                this.item[i].traceSphere(rayNow, rayHits, i);
                break;
            case BOX:
                this.item[i].traceBox(rayNow, rayHits, i);
                break;
            case CYLINDER:
                this.item[i].traceCylinder(rayNow, rayHits, i);
                break;
        }
    }
    rayHits.findNearestHit();
    if (rayHits.isShadowRay == true){
        return rayHits.iNearest;
    }
    var HitNow = rayHits.pierce[rayHits.iNearest];
    vec4.scale(HitNow.hitPt, rayNow.dir, HitNow.th);
    vec4.add(HitNow.hitPt, HitNow.hitPt, rayNow.orig);
    // Shading: find color at hit-point;
    return this.findShade(HitNow, rayHits.iNearest, rayNow, depth);
}
CScene.prototype.findShade = function(HitNow, hitNum, rayNow, depth){
    //if (depth <0){ return; }
    if (hitNum == 0){ return HitNow.colr; }
    var sumColr = vec3.create();
    var Astuff = this.item[HitNow.hitItem].matlNum;
    var reflectance = 1, transmittance = 0;
    
    var dir = vec4.scale(vec4.create(), HitNow.viewN, -1.0);
    var shadowP = vec4.scaleAndAdd(vec4.create(), HitNow.hitPt, dir, -2.92);
    this.item[HitNow.hitItem].findNormal(HitNow);       //surface normal
    
    
    if (Astuff >= 0 && light > 0){
        if (Astuff == 0 && HitNow.colr == this.item[HitNow.hitItem].gapColr){
            var MatlNow = this.matList[Astuff+1];
        }
        else if (HitNow.hitItem == 2 && Astuff == 3){
            var tot = Math.floor(HitNow.hitPt[0]) + Math.floor(HitNow.hitPt[1]) + Math.floor(HitNow.hitPt[2]);
            if(tot < 0) { var ans = -(tot % 2); }
            else { var ans = tot % 2; }
            if(ans) { var MatlNow = this.matList[5]; }
            else { var MatlNow = this.matList[Astuff]; }
        }
        else if (HitNow.hitItem == 3 && Astuff == 4){
            var tot = Math.floor(HitNow.hitPt[0]) + Math.floor(HitNow.hitPt[1]) + Math.floor(HitNow.hitPt[2]);
            if(tot < 0) { var ans = -(tot % 2); }
            else { var ans = tot % 2; }
            if(ans) { var MatlNow = this.matList[6]; }
            else { var MatlNow = this.matList[Astuff]; }
        }
        else{
            var MatlNow = this.matList[Astuff];
        }
        for (var i=0; i<this.lampCount; i++){
            if (this.lightList[i].enabled == true){
                var ambient = vec3.multiply(vec3.create(), this.lightList[i].I_a, MatlNow.kamb)
                vec3.add(sumColr, sumColr, ambient);
                
                if (i != 0 && shadowFlag > 0){      //shadow
                    var shadowDir = vec4.sub(vec4.create(), this.lightList[i].lightPos, shadowP);
                    vec4.normalize(shadowDir, shadowDir);
                    var shadowRay = new CRay();
                    vec4.copy(shadowRay.orig, shadowP);
                    vec4.copy(shadowRay.dir, shadowDir);
                    
                    var hitsList = new CHitList();
                    hitsList.initList(this.bkgndColr);
                    hitsList.initShadow();
                    var nearNum = this.trace(shadowRay, hitsList, depth);
                    
                    if (nearNum > 0){
                        vec3.scale(sumColr, sumColr, 0.3);
                        continue;
                    }
                }
                var L = vec4.sub(vec4.create(), this.lightList[i].lightPos, HitNow.hitPt);  //hitpt --> light
                vec4.normalize(L, L);
                var nDotL = Math.max(vec4.dot(L, HitNow.surfNorm), 0.0);
                var diffuse = vec3.multiply(vec3.create(), this.lightList[i].I_d, MatlNow.kdif);
                vec3.scaleAndAdd(sumColr, sumColr, diffuse, nDotL);
                
                var H = vec4.add(vec4.create(), L, HitNow.viewN);
                vec4.normalize(H, H);
                var nDotH = Math.max(vec4.dot(H, HitNow.surfNorm), 0.0);
                var e02 = nDotH * nDotH;
                var e04 = e02 * e02;
                var e08 = e04 * e04;
                var e16 = e08 * e08;
                var e32 = e16 * e16;
                var e64 = e32 * e32;
                var e128 = e64 * e64;
                var specular = vec3.multiply(vec3.create(), this.lightList[i].I_s, MatlNow.kspec);
                vec3.scaleAndAdd(sumColr, sumColr, specular, e128);
                
            }
        }
        
        if (depth > 0 && this.item[HitNow.hitItem].transparency > 0){        //refraction
            var cosI = vec4.dot(HitNow.surfNorm, dir);
            if (cosI > 0){
                var n1 = this.item[HitNow.hitItem].refractIndex;
                var n2 = 1.0;
                vec4.scale(HitNow.surfNorm, HitNow.surfNorm, -1.0);
            }
            else{
                var n1 = 1.0;
                var n2 = this.item[HitNow.hitItem].refractIndex;
                cosI = -cosI;
            }
            var ratio = n1 / n2;
            var sinT2 = ratio * ratio * (1.0 - cosI*cosI);
            var cosT2 = 1.0 - sinT2;
            if (cosT2 < 0){                 // Total internal reflection
                relectance = 1;
                transmittance = 0;
            }
            else{
                var cosT = Math.pow(cosT2, 0.5);
                var rn = (n1*cosI - n2*cosT) / (n1*cosI + n2*cosT);
                var rt = (n2*cosI - n1*cosT) / (n1*cosT + n2*cosI);
                rn *= rn;
                rt *= rt;
                reflectance = (rn + rt)*0.5;
                transmittance = 1 - reflectance;
                
                if (reflectance<=0 && transmittance<=0){
                    return sumColr;
                }
                if (transmittance > 0){
                    var tt = vec4.scale(vec4.create(), HitNow.surfNorm, (ratio*cosI-cosT));
                    var T = vec4.scaleAndAdd(vec4.create(), tt, dir, ratio);
                    vec4.normalize(T, T);
                    
                    var refractRay = new CRay();
                    var refractP = vec4.scaleAndAdd(vec4.create(), HitNow.hitPt, T, 0.001);
                    vec4.copy(refractRay.orig, refractP);
                    vec4.copy(refractRay.dir, T);
                    
                    var refractList = new CHitList();
                    refractList.initList(this.bkgndColr);
                    
                    var refractColr = this.trace(refractRay, refractList, depth-1);
                    vec3.scaleAndAdd(sumColr, sumColr, refractColr, this.item[HitNow.hitItem].transparency*transmittance);
                }
            }
        }
        
        if (depth > 0 && this.item[HitNow.hitItem].shininess > 0 && reflectance > 0){     //reflection
            var temp = vec4.dot(dir, HitNow.surfNorm)*2;
            var m = vec4.scale(vec4.create(), HitNow.surfNorm, temp);
            var r = vec4.sub(vec4.create(), dir, m);
            vec4.normalize(r, r);
            
            var reflectRay = new CRay();
            vec4.copy(reflectRay.orig, shadowP);
            vec4.copy(reflectRay.dir, r);
            
            var reflectList = new CHitList();
            reflectList.initList(this.bkgndColr);
            
            var reflectColr = this.trace(reflectRay, reflectList, depth-1);
            vec3.scaleAndAdd(sumColr, sumColr, reflectColr, this.item[HitNow.hitItem].shininess*reflectance);
            
        }
        
        return sumColr;
    }
    
        return HitNow.colr;

    
}
CScene.prototype.makeRayTracedImage = function(j,i){
    this.rayHits = new CHitList();
    this.rayHits.initList(this.bkgndColr);
    if (this.isAA == false){
        this.rayCam.uvnPix = vec4.fromValues(this.rayCam.iLeft+(i+0.5)*this.rayCam.ufrac, this.rayCam.iBot+(j+0.5)*this.rayCam.vfrac, this.rayCam.znear, 1.0);
        this.rayCam.makeEyeRay(this.rayNow);
        
        return this.trace(this.rayNow, this.rayHits, this.depthMax);
    }
    // Otherwise, use antialiasing!
    var uufrac = this.rayCam.ufrac / this.xSuperAA;
    var vvfrac = this.rayCam.vfrac / this.ySuperAA;
    var count = this.xSuperAA * this.ySuperAA;
    var countColr = vec3.create();
    for (var jj=0; jj<this.ySuperAA; jj++){
        for (var ii=0; ii<this.xSuperAA; ii++){
            this.rayHits = new CHitList();
            this.rayHits.initList(this.bkgndColr);
            this.rayCam.uvnPix = vec4.fromValues(this.rayCam.iLeft+i*this.rayCam.ufrac+(ii+Math.random())*uufrac, this.rayCam.iBot+j*this.rayCam.vfrac+(jj+Math.random())*vvfrac, this.rayCam.znear, 1.0);
            this.rayCam.makeEyeRay(this.rayNow);

            var colrTemp = this.trace(this.rayNow, this.rayHits, this.depthMax);
            vec3.add(countColr, countColr, colrTemp);
        }
    }
    vec3.scale(countColr, countColr, 1.0/count);
    return countColr;
}

//==============================================================================

function main() {
    //==============================================================================
    // Retrieve <canvas> element
    canvas = document.getElementById('webgl');
    
    browserResize();
    // Get the rendering context for WebGL
    gl = getWebGLContext(canvas);
    if (!gl) {
        console.log('Failed to get the rendering context for WebGL');
        return;
    }
    
    // Register the Mouse & Keyboard Event-handlers-------------------------------
    canvas.onmousedown	=	function(ev){myMouseDown( ev, gl, canvas) };
    canvas.onmousemove = 	function(ev){myMouseMove( ev, gl, canvas) };
    canvas.onmouseup = 		function(ev){myMouseUp(   ev, gl, canvas)};
    
    document.onkeydown = function(ev){ keydown(ev, gl); };
    window.addEventListener("keypress", myKeyPress, false);
    window.addEventListener("keydown", myKeyDown, false);
    window.addEventListener("keyup", myKeyUp, false);
    // END Mouse & Keyboard Event-Handlers-----------------------------------
    
    // Initialize shaders
    if (!initShaders(gl, VSHADER_SOURCE, FSHADER_SOURCE)) {
        console.log('Failed to intialize shaders.');
        return;
    }
    
    gl.depthFunc(gl.LESS);			 // WebGL default setting: (default)
    gl.enable(gl.DEPTH_TEST);
    
    // Set the vertex coordinates and color (the blue triangle is in the front)
    var n = initVertexBuffers(gl);
    if (n < 0) {
        console.log('Failed to specify the vertex infromation');
        return;
    }
    
    // Create, set uniform var to select fixed color vs texture map drawing:
    u_isTextureID = gl.getUniformLocation(gl.program, 'u_isTexture');
    u_flag = gl.getUniformLocation(gl.program, 'u_flag');
    u_isLighting = gl.getUniformLocation(gl.program, 'u_isLighting');
    u_Light0 = gl.getUniformLocation(gl.program, 'u_Light0');
    u_Light1 = gl.getUniformLocation(gl.program, 'u_Light1');
    if (!u_isTextureID || !u_flag || !u_isLighting || !u_Light0 || !u_Light1) {
        console.log('Failed to get the GPU storage location of uniform');
        return false;
    }
    
    // Specify the color for clearing <canvas>
    gl.clearColor(0.1, 0.1, 0.2, 1.0);
    
    // Get the storage locations of u_ViewMatrix and u_ProjMatrix variables
    u_eyePosWorld = gl.getUniformLocation(gl.program, 'u_eyePosWorld');
    u_NormalMatrix = gl.getUniformLocation(gl.program,'u_NormalMatrix');
    u_ViewMatrix = gl.getUniformLocation(gl.program, 'u_ViewMatrix');
    u_ProjMatrix = gl.getUniformLocation(gl.program, 'u_ProjMatrix');
    u_ModelMatrix = gl.getUniformLocation(gl.program, 'u_ModelMatrix');
    if (!u_ViewMatrix || !u_ProjMatrix || !u_ModelMatrix) {
        console.log('Failed to get u_ViewMatrix or u_ProjMatrix');
        return;
    }
    //  ... for Phong light source:
    u_Lamp0Pos  = gl.getUniformLocation(gl.program, 	'u_Lamp0Pos');
    u_Lamp0Amb  = gl.getUniformLocation(gl.program, 	'u_Lamp0Amb');
    u_Lamp0Diff = gl.getUniformLocation(gl.program, 	'u_Lamp0Diff');
    u_Lamp0Spec	= gl.getUniformLocation(gl.program,		'u_Lamp0Spec');
    if( !u_Lamp0Pos || !u_Lamp0Amb	|| !u_Lamp0Diff || !u_Lamp0Spec) {
        console.log('Failed to get the Lamp0 storage locations');
        return;
    }
    
    //  ... for Phong light source:
    u_Lamp1Pos  = gl.getUniformLocation(gl.program, 	'u_Lamp1Pos');
    u_Lamp1Amb  = gl.getUniformLocation(gl.program, 	'u_Lamp1Amb');
    u_Lamp1Diff = gl.getUniformLocation(gl.program, 	'u_Lamp1Diff');
    u_Lamp1Spec	= gl.getUniformLocation(gl.program,		'u_Lamp1Spec');
    if( !u_Lamp1Pos || !u_Lamp1Amb	|| !u_Lamp1Diff || !u_Lamp1Spec) {
        console.log('Failed to get the Lamp0 storage locations');
        return;
    }

    // ... for Phong material/reflectance:
    u_Ke = gl.getUniformLocation(gl.program, 'u_Ke');
    u_Ka = gl.getUniformLocation(gl.program, 'u_Ka');
    u_Kd = gl.getUniformLocation(gl.program, 'u_Kd');
    u_Ks = gl.getUniformLocation(gl.program, 'u_Ks');
    if(!u_Ke || !u_Ka || !u_Kd || !u_Ks) {
        console.log('Failed to get the Phong Reflectance storage locations');
        return;
    }
    
    // Position the first light source in World coords:
    gl.uniform3f(u_Lamp0Amb,  0.4, 0.4, 0.4);		// ambient
    gl.uniform3f(u_Lamp0Diff, 1.0, 1.0, 1.0);		// diffuse
    gl.uniform3f(u_Lamp0Spec, 1.0, 1.0, 1.0);		// Specular
    
    gl.uniform3f(u_Lamp1Amb,  0.05, 0.05, 0.05);		// ambient
    gl.uniform3f(u_Lamp1Diff, 0.5, 0.5, 0.5);		// diffuse
    gl.uniform3f(u_Lamp1Spec, 0.8, 0.8, 0.8);		// Specular
    
    projMatrix = new Matrix4();
    viewMatrix = new Matrix4();
    modelMatrix = new Matrix4();
    normalMatrix = new Matrix4();
    
    projMatrix.setPerspective(60, gl.drawingBufferWidth/2/gl.drawingBufferHeight, 1, 1000);
    gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
    
    createScene(); // use the function to create all the classes
    
    // Create, load, enable texture buffer object (TBO) in graphics hardware
    if (!initTextures(gl, n)) {
        console.log('Failed to intialize the texture object(s).');
        return;
    }
    
    //testQuaternions();		// test fcn at end of file
    
    draw(gl);
}

function createScene(){
    myCam = new CCamera(Math.PI / 3, gl.drawingBufferWidth/2/gl.drawingBufferHeight, 1, 1000);
    myCam.setFrustum();
    myScene = new CScene();
    
    myLight.lightPos = vec4.fromValues(g_EyeX, g_EyeY, g_EyeZ, 1.0);        //eye lamp
    myLight.I_a = vec3.fromValues(0.4, 0.4, 0.4);
    myLight.I_d = vec3.fromValues(1.0, 1.0, 1.0);
    myLight.I_s = vec3.fromValues(1.0, 1.0, 1.0);
    myScene.addLamp(myLight);
    HeadLight.I_a = vec3.fromValues(0.05, 0.05, 0.05);      //head lamp
    HeadLight.I_d = vec3.fromValues(0.5, 0.5, 0.5);
    HeadLight.I_s = vec3.fromValues(0.8, 0.8, 0.8);
    myScene.addLamp(HeadLight);
    
    myMatl.kemi = vec3.fromValues(0.0, 0.0, 0.0);          //grid line 0
    myMatl.kamb = vec3.fromValues(0.05375,  0.05,     0.06625);
    myMatl.kdif = vec3.fromValues(0.18275,  0.17,     0.22525);
    myMatl.kspec = vec3.fromValues(0.332741, 0.328634, 0.346435);
    myScene.addMatl(myMatl);
    myMatlB.kemi = vec3.fromValues(0.0, 0.0, 0.0);          //grid gap 1
    myMatlB.kamb = vec3.fromValues(2.5, 2.5, 2.5);
    myMatlB.kdif = vec3.fromValues(1.0, 1.0, 1.0);
    myMatlB.kspec = vec3.fromValues(1.0, 1.0, 1.0);
    myScene.addMatl(myMatlB);
    myMatlC.kemi = vec3.fromValues(0.0, 0.0, 0.0);          //sphere 2
    myMatlC.kamb = vec3.fromValues(0.1745,   0.01175,  0.01175);
    myMatlC.kdif = vec3.fromValues(1.0, 0.5, 0.5);
    myMatlC.kspec = vec3.fromValues(0.727811, 0.626959, 0.626959);
    myScene.addMatl(myMatlC);
    myMatlD.kemi = vec3.fromValues(0.0, 0.0, 0.0);           //sphereB 3
    myMatlD.kamb = vec3.fromValues(0.1, 0.18725, 0.1745);
    myMatlD.kdif = vec3.fromValues(0.396, 0.74151, 0.69102);
    myMatlD.kspec = vec3.fromValues(0.297254, 0.30829, 0.306678);
    myScene.addMatl(myMatlD);
    myMatlE.kemi = vec3.fromValues(0.0, 0.0, 0.0);          //box 4
    myMatlE.kamb = vec3.fromValues(0.2295,   0.08825,  0.0275);
    myMatlE.kdif = vec3.fromValues(0.9508,   0.9118,   0.566);
    myMatlE.kspec = vec3.fromValues(0.580594, 0.223257, 0.069570);
    myScene.addMatl(myMatlE);
    myMatlF.kemi = vec3.fromValues(0.0, 0.0, 0.0);          //sphere C 5
    myMatlF.kamb = vec3.fromValues(0.25,     0.20725,  0.20725);
    myMatlF.kdif = vec3.fromValues(1.0,      0.829,    0.829);
    myMatlF.kspec = vec3.fromValues(0.296648, 0.296648, 0.296648);
    myScene.addMatl(myMatlF);
    myMatlG.kemi = vec3.fromValues(0.0, 0.0, 0.0);          //sphere D 6
    myMatlG.kamb = vec3.fromValues(0.6, 0.1, 0.1);
    myMatlG.kdif = vec3.fromValues(0.98, 0.8, 0.68);
    myMatlG.kspec = vec3.fromValues(0.8, 0.8, 0.8);
    myScene.addMatl(myMatlG);
    
    myMatlH.kemi = vec3.fromValues(0.0, 0.0, 0.0);          //sphere E 7
    myMatlH.kamb = vec3.fromValues(0.05375,  0.05,     0.06625);
    myMatlH.kdif = vec3.fromValues(0.08275,  0.07,     0.02525);
    myMatlH.kspec = vec3.fromValues(0.332741, 0.328634, 0.346435);
    myScene.addMatl(myMatlH);
    
    myMatlI.kemi = vec3.fromValues(0.0, 0.0, 0.0);          // 8
    myMatlI.kamb = vec3.fromValues(0.329412, 0.223529, 0.027451);
    myMatlI.kdif = vec3.fromValues(0.780392, 0.568627, 0.113725);
    myMatlI.kspec = vec3.fromValues(0.992157, 0.941176, 0.807843);
    myScene.addMatl(myMatlI);
    myMatlJ.kemi = vec3.fromValues(0.0, 0.0, 0.0);          // 9
    myMatlJ.kamb = vec3.fromValues(0.329412, 0.223529, 0.027451);
    myMatlJ.kdif = vec3.fromValues(0.780392, 0.568627, 0.313725);
    myMatlJ.kspec = vec3.fromValues(0.992157, 0.941176, 0.807843);
    myScene.addMatl(myMatlJ);
    myMatlK.kemi = vec3.fromValues(0.0, 0.0, 0.0);          // 10
    myMatlK.kamb = vec3.fromValues(0.19125,  0.0735,   0.0225);
    myMatlK.kdif = vec3.fromValues(0.7038,   0.57048,  0.5828);
    myMatlK.kspec = vec3.fromValues(0.256777, 0.137622, 0.086014);
    myScene.addMatl(myMatlK);
    myMatlL.kemi = vec3.fromValues(0.0, 0.0, 0.0);          // 11
    myMatlL.kamb = vec3.fromValues(0.05375,  0.05,     0.06625);
    myMatlL.kdif = vec3.fromValues(0.58275,  0.57,     0.52525);
    myMatlL.kspec = vec3.fromValues(0.332741, 0.328634, 0.346435);
    myScene.addMatl(myMatlL);
    
    myMatlM.kemi = vec3.fromValues(0.0, 0.0, 0.0);          // 12
    myMatlM.kamb = vec3.fromValues(0.2295,   0.08825,  0.0275);
    myMatlM.kdif = vec3.fromValues(0.4508,   0.6118,   0.466);
    myMatlM.kspec = vec3.fromValues(0.580594, 0.223257, 0.0695701);
    myScene.addMatl(myMatlM);
    myMatlN.kemi = vec3.fromValues(0.0, 0.0, 0.0);          // 13
    myMatlN.kamb = vec3.fromValues(0.05375,  0.05,     0.06625);
    myMatlN.kdif = vec3.fromValues(0.58275,  0.57,     0.52525);
    myMatlN.kspec = vec3.fromValues(0.332741, 0.328634, 0.346435);
    myScene.addMatl(myMatlN);
    myMatlO.kemi = vec3.fromValues(0.0, 0.0, 0.0);          // 14
    myMatlO.kamb = vec3.fromValues(0.329412, 0.223529, 0.027451);
    myMatlO.kdif = vec3.fromValues(0.580392, 0.468627, 0.513725);
    myMatlO.kspec = vec3.fromValues(0.992157, 0.941176, 0.807843);
    myScene.addMatl(myMatlO);
    myMatlP.kemi = vec3.fromValues(0.0, 0.0, 0.0);          // 15
    myMatlP.kamb = vec3.fromValues(0.25,     0.20725,  0.20725);
    myMatlP.kdif = vec3.fromValues(1.0,      0.929,    0.929);
    myMatlP.kspec = vec3.fromValues(0.296648, 0.296648, 0.296648);
    myScene.addMatl(myMatlP);
    
    
    sceneChange();
    
}

function sceneChange(){
    // scene 1------------------------------------------------------------------------------------
    myScene.item = new Array();
    myGrid.matlNum = 0;
    myGrid.shininess = 0.0;
    myScene.addGeom(myGrid);
    myScene.itemCount = 1;
    if (scene > 0){
        mySphere.matlNum = 2;
        mySphere.shininess = 0.5;
        mySphere.transparency = 0;
        myScene.addGeom(mySphere);
        
        mySphereB.sphereColr = vec3.fromValues(0.95, 0.65, 0.88);
        mySphereB.matlNum = 3;
        mySphereB.shininess = 0;
        mySphereB.transparency = 0;
        myScene.addGeom(mySphereB);
        
        myBox.matlNum = 4;
        myBox.shininess = 0;
        myScene.addGeom(myBox);
        
        mySphereC.sphereColr = vec3.fromValues(1.0, 0.5, 0.5);
        mySphereC.matlNum = 5;
        mySphereC.shininess = 0.3;
        mySphereC.transparency = 0;
        myScene.addGeom(mySphereC);
        
        mySphereD.sphereColr = vec3.fromValues(1.0, 0.829, 0.829);
        mySphereD.matlNum = 6;
        mySphereD.shininess = 0.0;
        mySphereD.transparency = 0;
        myScene.addGeom(mySphereD);
        
        myCylinder.matlNum = 8;
        myCylinder.shininess = 0;
        myScene.addGeom(myCylinder);
        
        myCone.matlNum = 9;
        myCone.topRadius = 0.0;
        myScene.addGeom(myCone);
        
        myShade.matlNum = 10;
        myShade.topRadius = 0.5;
        myScene.addGeom(myShade);
        
        myPillar.matlNum = 11;
        myPillar.topRadius = 1.0;
        myScene.addGeom(myPillar);
        
        myBottom.matlNum = 7;
        myScene.addGeom(myBottom);
    }
    else{
        mySphereC.matlNum = 5;
        mySphereC.shininess = 0.5;
        mySphereC.transparency = 0;
        mySphereC.refractIndex = 2.0;
        myScene.addGeom(mySphereC);
        
        mySphere.matlNum = 7;
        mySphere.shininess = 0;
        mySphere.transparency = 0
        mySphere.refractIndex = 5;
        myScene.addGeom(mySphere);
        
        mySphereB.matlNum = 7;
        mySphereB.shininess = 0;
        mySphereB.transparency = 0
        mySphereB.refractIndex = 1.6;
        myScene.addGeom(mySphereB);
        
        myBottom.matlNum = 7;
        myScene.addGeom(myBottom);
        
        myBox.matlNum = 12;
        myScene.addGeom(myBox);
        
        myCube.matlNum = 13;
        myScene.addGeom(myCube);
        
        myCone.matlNum = 14;
        myCone.topRadius = 0.0;
        myScene.addGeom(myCone);
        
        mySphereD.matlNum = 7;
        mySphereD.shininess = 1.0;
        mySphereD.transparency = 1.0;
        mySphereD.refractIndex = 1.6;
        myScene.addGeom(mySphereD);
        
        mySphereE.matlNum = 15;
        mySphereE.shininess = 0.0;
        myScene.addGeom(mySphereE);
        
        mySphereF.matlNum = 15;
        myScene.addGeom(mySphereF);
        
        myCylinder.matlNum = 11;
        myScene.addGeom(myCylinder);
        
        mySphereG.matlNum = 3;
        myScene.addGeom(mySphereG);
        
        mySphereH.matlNum = 3;
        myScene.addGeom(mySphereH);
        
    }

}

function rayChange(){
    myRay.orig = vec4.fromValues(g_EyeX, g_EyeY, g_EyeZ, 1.0);
    myCam.setFrustum();
    myLight.lightPos = vec4.fromValues(g_EyeX, g_EyeY, g_EyeZ, 1.0);
}

function initVertexBuffers(gl) {
//==============================================================================
 
  myGrid.drawGrid();
    
  verticesTexCoords = new Float32Array([
    -1.0,  1.0, 0.0,  0.0, 0.0, 0.0,   	0.0, 1.0,   0.0, 0.0, 1.0,			// upper left corner,
    -1.0, -1.0, 0.0,  0.0, 0.0, 0.0,  	0.0, 0.0,	0.0, 0.0, 1.0,			// lower left corner,
     1.0,  1.0, 0.0,  0.0, 0.0, 0.0,  	1.0, 1.0,	0.0, 0.0, 1.0,			// upper right corner,
     1.0, -1.0, 0.0,  0.0, 0.0, 0.0,   	1.0, 0.0,	0.0, 0.0, 1.0,			// lower right corner.
  ]);
    
  ground = new Float32Array([
    -xymax, xymax, 0.01,   0.0, 0.0, 0.0,   -xymax, xymax,   0.0, 0.0, 1.0,
    -xymax,-xymax, 0.01,   0.0, 0.0, 0.0,   -xymax,-xymax,   0.0, 0.0, 1.0,
     xymax, xymax, 0.01,   0.0, 0.0, 0.0,    xymax, xymax,   0.0, 0.0, 1.0,
     xymax,-xymax, 0.01,   0.0, 0.0, 0.0,    xymax,-xymax,   0.0, 0.0, 1.0,
  ]);
    
  var sColr = vec3.fromValues(0.95, 0.65, 0.88);
  mySphere.drawSphere(sColr);
    
  myBox.drawBox();
    
  myCylinder.drawCylinder(1.0);
    
  mySiz = gndVerts.length + verticesTexCoords.length + sphVerts.length + sphVerts.length + cube.length + ground.length + sphVerts.length + sphVerts.length + cylVerts.length + cylVerts.length + cylVerts.length;
  var nn = mySiz / floatsPerVertex;

  var vertices = new Float32Array(mySiz);
  gndStart = 0;
  for(i=0,j=0; j< gndVerts.length; i++,j++) {
  	 vertices[i] = gndVerts[j];
  }
  TexStart = i;						// next we'll store the ground-plane;
  for(j=0; j< verticesTexCoords.length; i++, j++) {
     vertices[i] = verticesTexCoords[j];
  }
  sphStart = i;
  for(j=0; j< sphVerts.length; i++, j++) {
     vertices[i] = sphVerts[j];
  }
  var sColrB = vec3.fromValues(0.54, 0.53, 0.95);
  mySphereB.drawSphere(sColrB);
  sphStartB = i;
  for(j=0; j< sphVerts.length; i++, j++) {
     vertices[i] = sphVerts[j];
  }
  cubeStart = i;
  for(j=0; j< cube.length; i++, j++) {
     vertices[i] = cube[j];
  }
  groundStart = i;
  for(j=0; j< ground.length; i++, j++) {
     vertices[i] = ground[j];
  }
  var sColrC = vec3.fromValues(1.0, 0.5, 0.5);
  mySphereC.drawSphere(sColrC);
  sphStartC = i;
  for(j=0; j< sphVerts.length; i++, j++) {
     vertices[i] = sphVerts[j];
  }
  var sColrD = vec3.fromValues(1.0, 0.829, 0.829);
  mySphereD.drawSphere(sColrD);
  sphStartD = i;
  for(j=0; j< sphVerts.length; i++, j++) {
     vertices[i] = sphVerts[j];
  }
  cylStart = i;
  for(j=0; j< cylVerts.length; i++, j++) {
     vertices[i] = cylVerts[j];
  }
  myCone.drawCylinder(0.0);
  coneStart = i;
  for(j=0; j< cylVerts.length; i++, j++) {
     vertices[i] = cylVerts[j];
  }
  myShade.drawCylinder(0.5);
  shadeStart = i;
  for(j=0; j< cylVerts.length; i++, j++) {
     vertices[i] = cylVerts[j];
  }
  
  // Create a buffer object
  var vertexBufferID = gl.createBuffer();
  if (!vertexBufferID) {
    console.log('Failed to create the buffer object');
    return -1;
  }

  // Write vertex information to buffer object
  gl.bindBuffer(gl.ARRAY_BUFFER, vertexBufferID);
  gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.STATIC_DRAW);

  var FSIZE = vertices.BYTES_PER_ELEMENT;
  //---------------------------
  // Assign the buffer object to a_Position and enable the assignment
  var a_PositionID = gl.getAttribLocation(gl.program, 'a_Position');
  if(a_PositionID < 0) {
    console.log('Failed to get the storage location of a_Position');
    return -1;
  }
  gl.vertexAttribPointer(a_PositionID, 3, gl.FLOAT, false, FSIZE * 11, 0);
  gl.enableVertexAttribArray(a_PositionID);
  //---------------------------
  // Assign the buffer object to a_Color and enable the assignment
  var a_ColorID = gl.getAttribLocation(gl.program, 'a_Color');
  if(a_ColorID < 0) {
    console.log('Failed to get the storage location of a_Color');
    return -1;
  }
  gl.vertexAttribPointer(a_ColorID, 3, gl.FLOAT, false, FSIZE * 11, FSIZE * 3);
  gl.enableVertexAttribArray(a_ColorID);
  //---------------------------
  // Get the storage location of a_TexCoord
  var a_TexCoordID = gl.getAttribLocation(gl.program, 'a_TexCoord');
  if (a_TexCoordID < 0) {
    console.log('Failed to get the GPU storage location of a_TexCoord');
    return -1;
  }
  // Assign the buffer object to a_TexCoord variable
  gl.vertexAttribPointer(a_TexCoordID, 2, gl.FLOAT, false, FSIZE * 11, FSIZE * 6);
  gl.enableVertexAttribArray(a_TexCoordID);

  var a_Normal = gl.getAttribLocation(gl.program, 'a_Normal');
  if(a_Normal < 0) {
    console.log('Failed to get the storage location of a_Normal');
    return -1;
  }
  gl.vertexAttribPointer(a_Normal, 3, gl.FLOAT, false, FSIZE * 11, FSIZE * 8);
  gl.enableVertexAttribArray(a_Normal);
    
  return mySiz/floatsPerVertex;	// return # of vertices
}

function initTextures(gl, n) {
    //==============================================================================
    // set up the GPU to supply a texture image and pixel-by-pixel texture addresses
    // for our Fragment Shader.
    var textureID = gl.createTexture();   // Get GPU location for new texture map
    if (!textureID) {
        console.log('Failed to create the texture object on the GPU');
        return false;
    }
    
    // Get GPU location of a new uniform u_Sampler
    var u_SamplerID = gl.getUniformLocation(gl.program, 'u_Sampler');
    if (!u_SamplerID) {
        console.log('Failed to get the GPU location of u_Sampler');
        return false;
    }
    
    //myScene.rayNow.rayRotate(45, 1.0,0.0,0.0);
    //myGrid.rayRotate(45, 1.0,0.0,0.0);
    myPic.setTestPattern(1);				// fill it with an initial test-pattern.
//    myGrid.rayRotate(90, 1.0,0.0,0.0);
//    mat4.multiply(myRay.orig, myGrid.world2model, myRay.orig);
//    mat4.multiply(myRay.dir, myGrid.world2model, myRay.dir);
    //nuImag = myScene.makeRayTracedImage(nuImg);
    // 0 == colorful L-shaped pattern
    // 1 == uniform orange screen
    
    // Enable texture unit0 for our use
    gl.activeTexture(gl.TEXTURE0);
    // Bind our texture object (made at start of this fcn) to GPU's texture hdwe.
    gl.bindTexture(gl.TEXTURE_2D, textureID);
    // allocate memory and load the texture image into our texture object on GPU:
    gl.texImage2D(gl.TEXTURE_2D, 	//  'target'--the use of this texture
                  0, 							//  MIP-map level (default: 0)
                  gl.RGB, 				// GPU's data format (RGB? RGBA? etc)
                  myPic.xSiz,			// image width in pixels,
                  myPic.ySiz,			// image height in pixels,
                  0,							// byte offset to start of data
                  gl.RGB, 				// source/input data format (RGB? RGBA?)
                  gl.UNSIGNED_BYTE, 	// data type for each color channel
                  myPic.iBuf);	// data source.
    
    // Set the WebGL texture-filtering parameters
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    // Set the texture unit 0 to be driven by the sampler
    gl.uniform1i(u_SamplerID, 0);
    return true;									// done.
}

function refreshTextures(gl) {
    //==============================================================================
    // Modify/update the contents of the texture map(s) stored in the GPU;
    // copy current contents of CImgBuf object 'nuImg'  (see initTextures() above)
    // into the existing texture-map object stored in the GPU:
    
    gl.texSubImage2D(gl.TEXTURE_2D, 	//  'target'--the use of this texture
                     0, 							//  MIP-map level (default: 0)
                     0,0,						// xoffset, yoffset (shifts the image)
                     myPic.xSiz,			// image width in pixels,
                     myPic.ySiz,			// image height in pixels,
                     gl.RGB, 				// source/input data format (RGB? RGBA?)
                     gl.UNSIGNED_BYTE, 	// data type for each color channel
                     myPic.iBuf);	// data source.
}

function draw(gl) {
//==============================================================================
  
  // Clear <canvas> color AND DEPTH buffer
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

//------------------------------------------
// Draw in the LEFT viewport
//------------------------------------------
	// CHANGE from our default viewport:
	// gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight);
	// to a smaller one:
  gl.viewport(0,  														// Viewport lower-left corner
							0,										   // (x,y) location(in pixels)
  						gl.drawingBufferWidth/2, 				// viewport width, height.
  						gl.drawingBufferHeight);
  // Set the matrix to be used for to set the camera view
  viewMatrix.setLookAt(g_EyeX, g_EyeY, g_EyeZ, g_AtX, g_AtY, g_AtZ, 0, 1, 0);
  gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
    
  gl.uniform4f(u_eyePosWorld, g_EyeX, g_EyeY, g_EyeZ, 1);
  gl.uniform4f(u_Lamp0Pos,  g_EyeX, g_EyeY, g_EyeZ, 1.0);
  gl.uniform4f(u_Lamp1Pos,  lampX, lampY, lampZ, 1.0);
  HeadLight.lightPos = vec4.fromValues(lampX, lampY, lampZ, 1.0);
  // select fixed-color drawing:
  gl.uniform1i(u_flag, 0);
  gl.uniform1i(u_isTextureID, 0);
  gl.uniform1i(u_isLighting, light);
  gl.uniform1i(u_Light0, light0);
  gl.uniform1i(u_Light1, light1);
  drawMyScene(gl);
    
    
//------------------------------------------
// Draw in the RIGHT viewport:
//------------------------------------------
  gl.viewport(gl.drawingBufferWidth/2, 				// Viewport lower-left corner
							0, 															// location(in pixels)
  						gl.drawingBufferWidth/2, 				// viewport width, height.
  						gl.drawingBufferHeight);
  gl.uniform1i(u_flag, 1);
  gl.uniform1i(u_isTextureID, 1);
  gl.drawArrays(gl.TRIANGLE_STRIP,							// use this drawing primitive, and
                  TexStart/floatsPerVertex,	// start at this vertex number, and
                  verticesTexCoords.length/floatsPerVertex);		// draw this many vertices
}

function drawMyScene(myGL) {
//===============================================================================
 // Rotate to make a new set of 'world' drawing axes: 
 // old one had "+y points upwards", but
  modelMatrix.setTranslate(0,0,0);
  quatMatrix.setFromQuat(qTot.x, qTot.y, qTot.z, qTot.w);	// Quaternion-->Matrix
  modelMatrix.concat(quatMatrix);	// apply that matrix.
  //Grid------------------------------------------
    // Set the Phong materials' reflectance:
  gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
  gl.uniform3f(u_Ka, 0.05375,  0.05,     0.06625);				// Ka ambient
  gl.uniform3f(u_Kd, 0.18275,  0.17,     0.22525);				// Kd diffuse
  gl.uniform3f(u_Ks, 0.332741, 0.328634, 0.346435);				// Ks specular
    
  pushMatrix(modelMatrix);
  modelMatrix.translate(0.0, -1.0, myGrid.zVal);
  modelMatrix.rotate(90, 1,0,0);	// new one has "+z points upwards",
  modelMatrix.scale(2.0, 2.0, 2.0);
    
  // Calculate the matrix to transform the normal based on the model matrix
  normalMatrix.setInverseOf(modelMatrix);
  normalMatrix.transpose();
  gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
    
  myGrid.model2world = mat4.clone(modelMatrix.elements);
  mat4.invert(myGrid.world2model, myGrid.model2world);
  myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
  // Now, using these drawing axes, draw our ground plane: 
  myGL.drawArrays(myGL.LINES,							// use this drawing primitive, and
  							gndStart/floatsPerVertex,	// start at this vertex number, and
  							gndVerts.length/floatsPerVertex);		// draw this many vertices
    
  gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
  gl.uniform3f(u_Ka, 2.5, 2.5, 2.5);				// Ka ambient
  gl.uniform3f(u_Kd, 1.0, 1.0, 1.0);				// Kd diffuse
  gl.uniform3f(u_Ks, 1.0, 1.0, 1.0);				// Ks specular
  myGL.drawArrays(myGL.TRIANGLE_STRIP, groundStart/floatsPerVertex, ground.length/floatsPerVertex);
// scene 1------------------------------------------------------------------------------------
  if (scene > 0){
      //Sphere1------------------------------------------
      // Set the Phong materials' reflectance:
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.1745,   0.01175,  0.01175);				// Ka ambient
      gl.uniform3f(u_Kd, 1.0, 0.5, 0.5);				            // Kd diffuse
      gl.uniform3f(u_Ks, 0.727811, 0.626959, 0.626959);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(2.7, 0.5, -6.4);
      modelMatrix.scale(1.5, 1.5, 1.5);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      mySphere.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(mySphere.world2model, mySphere.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      sphStart/floatsPerVertex,	// start at this vertex number, and
                      sphVerts.length/floatsPerVertex);
      //Sphere2------------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.1,      0.18725,  0.1745);				// Ka ambient
      gl.uniform3f(u_Kd, 0.396,    0.74151,  0.69102);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.297254, 0.30829,  0.306678);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(-2.5, 1.0, -11.5);
      modelMatrix.scale(2.2, 2.2, 2.2);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      mySphereB.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(mySphereB.world2model, mySphereB.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      sphStartB/floatsPerVertex,	// start at this vertex number, and
                      sphVerts.length/floatsPerVertex);
      //Sphere3------------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.25,     0.20725,  0.20725);				// Ka ambient
      gl.uniform3f(u_Kd, 1.0,      0.829,    0.829);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.296648, 0.296648, 0.296648);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(0.63, 3.2, -10.5);
      modelMatrix.rotate(60, 0,0,1);
      modelMatrix.scale(1.0, 2.3, 1.0);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      mySphereC.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(mySphereC.world2model, mySphereC.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      sphStartC/floatsPerVertex,	// start at this vertex number, and
                      sphVerts.length/floatsPerVertex);
      //Sphere4------------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.6, 0.1, 0.1);				// Ka ambient
      gl.uniform3f(u_Kd, 0.98, 0.8, 0.68);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.8, 0.8, 0.8);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(0.1, 0.0, -4.0);
      modelMatrix.scale(0.35, 0.35, 0.35);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      mySphereD.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(mySphereD.world2model, mySphereD.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      sphStartD/floatsPerVertex,	// start at this vertex number, and
                      sphVerts.length/floatsPerVertex);
      
      //Box------------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive, MATL_COPPER_SHINY
      gl.uniform3f(u_Ka, 0.2295,   0.08825,  0.0275);				// Ka ambient
      gl.uniform3f(u_Kd, 0.9508,   0.9118,   0.566);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.580594, 0.223257, 0.0695701);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(1.2, 0.6, -10.7);
      modelMatrix.rotate(10, 0,1,0);
      modelMatrix.scale(2.5, 2.5, 2.5);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      myBox.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(myBox.world2model, myBox.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLES, 				// use this drawing primitive, and
                      cubeStart/floatsPerVertex,	// start at this vertex number, and
                      cube.length/floatsPerVertex);
      
      //cylinder------------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive, MATL_COPPER_SHINY
      gl.uniform3f(u_Ka, 0.329412, 0.223529, 0.027451);				// Ka ambient
      gl.uniform3f(u_Kd, 0.780392, 0.568627, 0.113725);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.992157, 0.941176, 0.807843);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(0.1, -0.3, -4.1);
      modelMatrix.rotate(90, 1,0,0);
      modelMatrix.scale(0.9, 0.9, 0.5);
      
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      myCylinder.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(myCylinder.world2model, myCylinder.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      cylStart/floatsPerVertex,	// start at this vertex number, and
                      cylVerts.length/floatsPerVertex);
      
      //cone------------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive, MATL_COPPER_SHINY
      gl.uniform3f(u_Ka, 0.329412, 0.223529, 0.027451);				// Ka ambient
      gl.uniform3f(u_Kd, 0.780392, 0.568627, 0.313725);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.992157, 0.941176, 0.807843);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(-2.44, 2.8, -11.5);
      modelMatrix.rotate(-90, 1,0,0);
      modelMatrix.scale(1.6, 1.6, 2.8);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      myCone.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(myCone.world2model, myCone.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      coneStart/floatsPerVertex,	// start at this vertex number, and
                      cylVerts.length/floatsPerVertex);
      
      //shade------------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive, MATL_COPPER_SHINY
      gl.uniform3f(u_Ka, 0.19125,  0.0735,   0.0225);				// Ka ambient
      gl.uniform3f(u_Kd, 0.7038,   0.57048,  0.5828);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.256777, 0.137622, 0.086014);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(-2.3, 0.5, -5.2);
      modelMatrix.rotate(-90, 1,0,0);
      modelMatrix.scale(1.0, 1.0, 1.1);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      myShade.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(myShade.world2model, myShade.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      shadeStart/floatsPerVertex,	// start at this vertex number, and
                      cylVerts.length/floatsPerVertex);
      
      //pillar------------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.05375,  0.05,     0.06625);				// Ka ambient
      gl.uniform3f(u_Kd, 0.58275,  0.57,     0.52525);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.332741, 0.328634, 0.346435);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(-2.3, -0.5, -5.2);
      modelMatrix.rotate(-90, 1,0,0);
      modelMatrix.scale(0.3, 0.3, 1.0);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      myPillar.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(myPillar.world2model, myPillar.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      cylStart/floatsPerVertex,	// start at this vertex number, and
                      cylVerts.length/floatsPerVertex);
      
      //bottom------------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.05375,  0.05,     0.06625);				// Ka ambient
      gl.uniform3f(u_Kd, 0.08275,  0.07,     0.02525);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.332741, 0.328634, 0.346435);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(-2.3, -0.5, -5.2);
      modelMatrix.scale(0.9, 0.3, 0.9);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      myBottom.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(myBottom.world2model, myBottom.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLES, 				// use this drawing primitive, and
                      cubeStart/floatsPerVertex,	// start at this vertex number, and
                      cube.length/floatsPerVertex);
      
    }
    
// scene 2------------------------------------------------------------------------------------
  else{
      //Sphere1------------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.25,     0.20725,  0.20725);				// Ka ambient
      gl.uniform3f(u_Kd, 1.0,      0.829,    0.829);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.296648, 0.296648, 0.296648);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(-1.5, 0.2, 0.9);
      modelMatrix.rotate(15, 0,0,1);
      pushMatrix(modelMatrix);
      modelMatrix.translate(0.0, 1.0, -6.5);
      modelMatrix.scale(1.3, 1.3, 1.3);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      mySphereC.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(mySphereC.world2model, mySphereC.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      sphStartC/floatsPerVertex,	// start at this vertex number, and
                      sphVerts.length/floatsPerVertex);
      //Sphere2------------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.05375,  0.05,     0.06625);				// Ka ambient
      gl.uniform3f(u_Kd, 0.08275,  0.07,     0.02525);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.332741, 0.328634, 0.346435);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(-1.1, 2.0, -6.2);
      modelMatrix.scale(0.5, 0.5, 0.5);

      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      mySphere.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(mySphere.world2model, mySphere.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      sphStart/floatsPerVertex,	// start at this vertex number, and
                      sphVerts.length/floatsPerVertex);
      //Sphere3------------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.05375,  0.05,     0.06625);				// Ka ambient
      gl.uniform3f(u_Kd, 0.08275,  0.07,     0.02525);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.332741, 0.328634, 0.346435);				// Ks specular
      
      modelMatrix = popMatrix();
      modelMatrix.translate(1.14, 2.0, -6.2);
      modelMatrix.scale(0.5, 0.5, 0.5);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      mySphereB.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(mySphereB.world2model, mySphereB.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      sphStartB/floatsPerVertex,	// start at this vertex number, and
                      sphVerts.length/floatsPerVertex);
      
      //Sphere4-----------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.05375,  0.05,     0.06625);				// Ka ambient
      gl.uniform3f(u_Kd, 0.08275,  0.07,     0.02525);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.332741, 0.328634, 0.346435);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(0.0, 0.0, 0.0);
      modelMatrix.rotate(20, 0,0,1);
      pushMatrix(modelMatrix);
      modelMatrix.translate(0.3, 1.4, -1.6);
      modelMatrix.scale(0.5, 0.5, 0.5);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      mySphereD.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(mySphereD.world2model, mySphereD.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      sphStartD/floatsPerVertex,	// start at this vertex number, and
                      sphVerts.length/floatsPerVertex);
      //Sphere5-----------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.25,     0.20725,  0.20725);				// Ka ambient
      gl.uniform3f(u_Kd, 1.0,      0.929,    0.929);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.296648, 0.296648, 0.296648);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(-0.2, 1.7, -1.5);
      modelMatrix.scale(0.19, 0.19, 0.19);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      mySphereE.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(mySphereE.world2model, mySphereE.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      sphStartC/floatsPerVertex,	// start at this vertex number, and
                      sphVerts.length/floatsPerVertex);
      //Sphere6-----------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.25,     0.20725,  0.20725);				// Ka ambient
      gl.uniform3f(u_Kd, 1.0,      0.929,    0.929);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.296648, 0.296648, 0.296648);				// Ks specular
      
      modelMatrix = popMatrix();
      modelMatrix.translate(0.8, 1.7, -1.6);
      modelMatrix.scale(0.19, 0.19, 0.19);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      mySphereF.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(mySphereF.world2model, mySphereF.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      sphStartC/floatsPerVertex,	// start at this vertex number, and
                      sphVerts.length/floatsPerVertex);
      
      //Bottom-----------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.05375,  0.05,     0.06625);				// Ka ambient
      gl.uniform3f(u_Kd, 0.08275,  0.07,     0.02525);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.332741, 0.328634, 0.346435);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(4.0, -0.8, -12.2);
      modelMatrix.scale(3.4, 0.3, 3.4);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      myBottom.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(myBottom.world2model, myBottom.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLES, 				// use this drawing primitive, and
                      cubeStart/floatsPerVertex,	// start at this vertex number, and
                      cube.length/floatsPerVertex);
      //Pillar-----------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive, MATL_COPPER_SHINY
      gl.uniform3f(u_Ka, 0.2295,   0.08825,  0.0275);				// Ka ambient
      gl.uniform3f(u_Kd, 0.4508,   0.6118,   0.466);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.580594, 0.223257, 0.0695701);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(4.0, 1.0, -12.7);
      modelMatrix.rotate(0, 0,1,0);
      modelMatrix.scale(2.0, 3.4, 2.0);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      myBox.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(myBox.world2model, myBox.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLES, 				// use this drawing primitive, and
                      cubeStart/floatsPerVertex,	// start at this vertex number, and
                      cube.length/floatsPerVertex);
      //Support-----------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.05375,  0.05,     0.06625);				// Ka ambient
      gl.uniform3f(u_Kd, 0.58275,  0.57,     0.52525);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.332741, 0.328634, 0.346435);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(3.8, 2.8, -12.2);
      modelMatrix.rotate(10, 0,1,0);
      modelMatrix.scale(2.2, 0.3, 2.2);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      myCube.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(myCube.world2model, myCube.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLES, 				// use this drawing primitive, and
                      cubeStart/floatsPerVertex,	// start at this vertex number, and
                      cube.length/floatsPerVertex);
      //Roof-----------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive, MATL_COPPER_SHINY
      gl.uniform3f(u_Ka, 0.329412, 0.223529, 0.027451);				// Ka ambient
      gl.uniform3f(u_Kd, 0.580392, 0.468627, 0.513725);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.992157, 0.941176, 0.807843);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(4.0, 2.86, -12.5);
      modelMatrix.rotate(-90, 1,0,0);
      modelMatrix.scale(1.2, 1.2, 3.3);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      myCone.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(myCone.world2model, myCone.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      coneStart/floatsPerVertex,	// start at this vertex number, and
                      cylVerts.length/floatsPerVertex);
      
      //cylinder------------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.05375,  0.05,     0.06625);				// Ka ambient
      gl.uniform3f(u_Kd, 0.58275,  0.57,     0.52525);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.332741, 0.328634, 0.346435);				// Ks specular
         
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(2.6, -0.6, -6.3);
      modelMatrix.rotate(90, 1,0,0);
      modelMatrix.scale(1.2, 1.0, 0.3);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      myCylinder.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(myCylinder.world2model, myCylinder.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      cylStart/floatsPerVertex,	// start at this vertex number, and
                      cylVerts.length/floatsPerVertex);
      //sphere7------------------------------------------
      gl.uniform3f(u_Ke, 0.0, 0.0, 0.0);				// Ke emissive
      gl.uniform3f(u_Ka, 0.1,      0.18725,  0.1745);				// Ka ambient
      gl.uniform3f(u_Kd, 0.396,    0.74151,  0.69102);				// Kd diffuse
      gl.uniform3f(u_Ks, 0.297254, 0.30829,  0.306678);				// Ks specular
      
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(2.5, -0.2, -6.3);
      modelMatrix.scale(0.62, 0.5, 0.5);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      mySphereG.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(mySphereG.world2model, mySphereG.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      sphStart/floatsPerVertex,	// start at this vertex number, and
                      sphVerts.length/floatsPerVertex);
      //sphere8------------------------------------------    
      modelMatrix = popMatrix();
      pushMatrix(modelMatrix);
      modelMatrix.translate(3.2, -0.2, -6.2);
      modelMatrix.scale(0.3, 0.3, 0.3);
      
      normalMatrix.setInverseOf(modelMatrix);
      normalMatrix.transpose();
      gl.uniformMatrix4fv(u_NormalMatrix, false, normalMatrix.elements);
      
      mySphereH.model2world = mat4.clone(modelMatrix.elements);
      mat4.invert(mySphereH.world2model, mySphereH.model2world);
      myGL.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
      myGL.drawArrays(myGL.TRIANGLE_STRIP, 				// use this drawing primitive, and
                      sphStart/floatsPerVertex,	// start at this vertex number, and
                      sphVerts.length/floatsPerVertex);
      
  }
    
}


function CImgBuf(wide, tall) {
    //==============================================================================
    // Construct an 'image-buffer' object to hold a floating-point ray-traced image.
    //  Contains BOTH
    //	iBuf -- 2D array of 8-bit RGB pixel values we can display on-screen, AND
    //	fBuf -- 2D array of floating-point RGB pixel values we CAN'T display.
    //			--Both buffers hold the same numbers of pixel values (xSiz,ySiz,pixSiz)
    //			--imgBuf.int2float() copies/converts current iBuf contents to fBuf
    //			--imgBuf.float2int() copies/converts current fBuf contents to iBuf
    //	WHY?
    //	--Our ray-tracer computes floating-point light amounts(e.g. radiance L) //    but neither our display nor our WebGL texture-map buffers can accept
    //		images with floating-point pixel values.
    //	--You will NEED all those floating-point values for many applications!
    // Stay simple in early versions of your ray-tracer: keep 0.0 <= RGB < 1.0,
    // but later you can modify your ray-tracer
    // to use radiometric units of Radiance (watts/(steradians*meter^2), or convert
    // to use photometric units of luminance (lumens/(steradians*meter^2 aka cd/m^2) // to compute in physically verifiable units of visible light.
    
    this.xSiz = wide;							// image width in pixels
    this.ySiz =	tall;							// image height in pixels
    this.pixSiz = 3;							// pixel size (3 for RGB, 4 for RGBA, etc)
    this.iBuf = new Uint8Array(  this.xSiz * this.ySiz * this.pixSiz);
    this.fBuf = new Float32Array(this.xSiz * this.ySiz * this.pixSiz);
}

CImgBuf.prototype.setTestPattern = function(pattNum) {
    //==============================================================================
    // Replace current 8-bit RGB contents of 'imgBuf' with a colorful pattern
    // 2D color image:  8-bit unsigned integers in a 256*256*3 array
    // to store r,g,b,r,g,b integers (8-bit)
    // In WebGL texture map sizes MUST be a power-of-two (2,4,8,16,32,64,...4096)
    // with origin at lower-left corner
    // (declare some global vars:)
    
    // use local vars to set the array's contents.
    for(var j=0; j< this.ySiz; j++) {						// for the j-th row of pixels
        for(var i=0; i< this.xSiz; i++) {					// and the i-th pixel on that row,
            var idx = (j*this.xSiz + i)*this.pixSiz;// Array index at pixel (i,j)
            switch(pattNum) {
                case 0:
                    var colr = myScene.makeRayTracedImage(j,i);
                    this.fBuf[idx   ] = colr[0];
                    this.fBuf[idx +1] = colr[1];
                    this.fBuf[idx +2] = colr[2];
                    break;
                case 1:
                    this.fBuf[idx   ] = 1.0;	// bright orange
                    this.fBuf[idx +1] = 1.0;
                    this.fBuf[idx +2] = 0.5;
                    break;
                default:
                    console.log("imgBuf.setTestPattern() says: WHUT!?");
                    break;
            }
        }
    }
    this.float2int();		// fill the floating-point buffer with same test pattern.
}

CImgBuf.prototype.int2float = function() {
    //==============================================================================
    // Convert current integerRGB image in iBuf into floating-point RGB image in fBuf
    for(var j=0; j< this.ySiz; j++) {		// for each scanline
        for(var i=0; i< this.xSiz; i++) {		// for each pixel on that scanline
            var idx = (j*this.xSiz + i)*this.pixSiz;// Find array index at pixel (i,j)
            // convert integer 0 <= RGB <= 255 to floating point 0.0 <= R,G,B <= 1.0
            this.fBuf[idx   ] = this.iBuf[idx   ] / 255.0;	// red
            this.fBuf[idx +1] = this.iBuf[idx +1] / 255.0;	// grn
            this.fBuf[idx +2] = this.iBuf[idx +2] / 255.0;	// blu
            
        }
    }
}

CImgBuf.prototype.float2int = function() {
    //==============================================================================
    // Convert current floating-point RGB image in fBuf into integerRGB image in iBuf
    for(var j=0; j< this.ySiz; j++) {		// for each scanline
        for(var i=0; i< this.xSiz; i++) {	// for each pixel on that scanline
            var idx = (j*this.xSiz + i)*this.pixSiz;// Find array index at pixel (i,j)
            // find 'clamped' color values that stay >=0.0 and <=1.0:
            var rval = Math.min(1.0, Math.max(0.0, this.fBuf[idx   ]))
            var gval = Math.min(1.0, Math.max(0.0, this.fBuf[idx +1]));
            var bval = Math.min(1.0, Math.max(0.0, this.fBuf[idx +2]));
            // Divide [0,1] span into 256 equal-sized parts: e.g.  Math.floor(rval*256)
            // In the rare case when rval==1.0 you get unwanted '256' result that won't
            // fit into the 8-bit RGB values.  Fix it with Math.max():
            this.iBuf[idx   ] = Math.min(255,Math.floor(rval*256.0));	// red
            this.iBuf[idx +1] = Math.min(255,Math.floor(gval*256.0));	// grn
            this.iBuf[idx +2] = Math.min(255,Math.floor(bval*256.0));	// blu
            
        }
    }
}



//===================Mouse and Keyboard event-handling Callbacks================
//==============================================================================
function myMouseDown(ev, gl, canvas) {
    //==============================================================================
    // Called when user PRESSES down any mouse button;
    // 									(Which button?    console.log('ev.button='+ev.button);   )
    // 		ev.clientX, ev.clientY == mouse pointer location, but measured in webpage
    //		pixels: left-handed coords; UPPER left origin; Y increases DOWNWARDS (!)
    
    // Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
    var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
    var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
    var yp = canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
    //  console.log('myMouseDown(pixel coords): xp,yp=\t',xp,',\t',yp);
    
    // Convert to Canonical View Volume (CVV) coordinates too:
    // MODIFIED for side-by-side display: find position within the LEFT-side CVV
    var x = (xp - canvas.width/4)  / 		// move origin to center of LEFT viewport,
    (canvas.width/4);			// normalize canvas to -1 <= x < +1,
    var y = (yp - canvas.height/2) /		//										 -1 <= y < +1.
    (canvas.height/2);
    //	console.log('myMouseDown(CVV coords  ):  x, y=\t',x,',\t',y);
    
    isDrag = true;											// set our mouse-dragging flag
    xMclik = x;													// record where mouse-dragging began
    yMclik = y;
};

function myMouseMove(ev, gl, canvas) {
    //==============================================================================
    // Called when user MOVES the mouse with a button already pressed down.
    // 									(Which button?   console.log('ev.button='+ev.button);    )
    // 		ev.clientX, ev.clientY == mouse pointer location, but measured in webpage
    //		pixels: left-handed coords; UPPER left origin; Y increases DOWNWARDS (!)
    
    if(isDrag==false) return;				// IGNORE all mouse-moves except 'dragging'
    
    // Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
    var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
    var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
    var yp = canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
    //  console.log('myMouseMove(pixel coords): xp,yp=\t',xp,',\t',yp);
    
    // Convert to Canonical View Volume (CVV) coordinates too:
    // MODIFIED for side-by-side display: find position within the LEFT-side CVV
    var x = (xp - canvas.width/4)  / 		// move origin to center of LEFT viewport,
    (canvas.width/4);			// normalize canvas to -1 <= x < +1,
    var y = (yp - canvas.height/2) /		//										 -1 <= y < +1.
    (canvas.height/2);
    //	console.log('myMouseMove(CVV coords  ):  x, y=\t',x,',\t',y);
    
    // find how far we dragged the mouse:
    xMdragTot += (x - xMclik);					// Accumulate change-in-mouse-position,&
    yMdragTot += (y - yMclik);
    
    // AND use any mouse-dragging we found to update quaternions qNew and qTot.
    dragQuat(x - xMclik, y - yMclik);
    draw(gl);
    
    xMclik = x;													// Make next drag-measurement from here.
    yMclik = y;
    
//    var dist = Math.sqrt(xMdragTot*xMdragTot + yMdragTot*yMdragTot);
//    myGrid.rayRotate(dist*120.0, -yMdragTot+0.0001, xMdragTot+0.0001, 0.0);
//    mat4.multiply(myRay.orig, myGrid.world2model, myRay.orig);
//    //console.log('ml', myGrid.world2model);
//    //console.log('or', myRay.orig);
//    draw(gl);
};

function myMouseUp(ev, gl, canvas) {
    //==============================================================================
    // Called when user RELEASES mouse button pressed previously.
    // 									(Which button?   console.log('ev.button='+ev.button);    )
    // 		ev.clientX, ev.clientY == mouse pointer location, but measured in webpage
    //		pixels: left-handed coords; UPPER left origin; Y increases DOWNWARDS (!)
    
    // Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
    var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
    var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
    var yp = canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
    //  console.log('myMouseUp  (pixel coords): xp,yp=\t',xp,',\t',yp);
    
    // Convert to Canonical View Volume (CVV) coordinates too:
    // MODIFIED for side-by-side display: find position within the LEFT-side CVV
    var x = (xp - canvas.width/4)  / 		// move origin to center of LEFT viewport,
    (canvas.width/4);			// normalize canvas to -1 <= x < +1,
    var y = (yp - canvas.height/2) /		//										 -1 <= y < +1.
    (canvas.height/2);
    //console.log('myMouseUp  (CVV coords  ):  x, y=\t',x,',\t',y);
    
    isDrag = false;											// CLEAR our mouse-dragging flag, and
    // accumulate any final bit of mouse-dragging we did:
    xMdragTot += (x - xMclik);
    yMdragTot += (y - yMclik);
    
    // AND use any mouse-dragging we found to update quaternions qNew and qTot;
    dragQuat(x - xMclik, y - yMclik);
    draw(gl);
    //console.log('myMouseUp: xMdragTot,yMdragTot =',xMdragTot,',\t',yMdragTot);
//    var dist = Math.sqrt(xMdragTot*xMdragTot + yMdragTot*yMdragTot);
//    myGrid.rayRotate(dist*120.0, -yMdragTot+0.0001, xMdragTot+0.0001, 0.0);
//    mat4.multiply(myRay.orig, myGrid.world2model, myRay.orig);
//    draw(gl);
};

function myKeyDown(ev) {
    //==============================================================================
//    switch(ev.keyCode){
//        case 187:       //+
//            zVal += 0.1;
//            break;
//        case 189:       //-
//            zVal -= 0.1;
//            break;
//    }
}

function myKeyUp(ev) {
    //==============================================================================
    // Called when user releases ANY key on the keyboard; captures scancodes well
    // You probably don't want to use this ('myKeyDown()' explains why); you'll find
    // myKeyPress() can handle nearly all your keyboard-interface needs.
    /*
     console.log('myKeyUp()--keyCode='+ev.keyCode+' released.');
     */
}

function myKeyPress(ev) {
    //===============================================================================
    // Best for capturing alphanumeric keys and key-combinations such as
    // CTRL-C, alt-F, SHIFT-4, etc.  Use this instead of myKeyDown(), myKeyUp() if
    // you don't need to respond separately to key-down and key-up events.
    
    /*
     // Report EVERYTHING about this pressed key in the console:
     console.log('myKeyPress():keyCode='+ev.keyCode  +', charCode=' +ev.charCode+
     ', shift='    +ev.shiftKey + ', ctrl='    +ev.ctrlKey +
     ', altKey='   +ev.altKey   +
     ', metaKey(Command key or Windows key)='+ev.metaKey);
     */
    myChar = String.fromCharCode(ev.keyCode);	//	convert code to character-string
    
    var myCanvas = document.getElementById('webgl');	// get current canvas
    var myGL = getWebGLContext(myCanvas);				// and its current context:
    
    switch(myChar){
        case 't':
        case 'T':
            if(isTracing == 0) isTracing = 1;
            else isTracing = 0;
            myPic.setTestPattern(isTracing);
            
            refreshTextures(myGL);
            draw(myGL)
            break;
        case 'p':
            if(myScene.isAA == true) myScene.isAA = false;
            else myScene.isAA = true;
            break;
        case 'o':
            if(light == 1) light = 0;
            else light = 1;
            draw(myGL)
            break;
        case 'k':
            if(light0 == 1) {light0 = 0; myLight.enabled = false;}
            else {light0 = 1; myLight.enabled = true;}
            draw(myGL)
            break;
        case 'l':
            if(light1 == 1) {light1 = 0; HeadLight.enabled = false;}
            else {light1 = 1; HeadLight.enabled = true;}
            draw(myGL)
            break;
        case 'm':
            if(shadowFlag == 1) shadowFlag = 0;
            else shadowFlag = 1;
            break;
        case 'n':
            if(scene == 1) scene = 0;
            else scene = 1;
            draw(myGL);
            sceneChange();
            break;
       
    }
    

}

function keydown(ev, gl) {
    //------------------------------------------------------
    //HTML calls this'Event handler' or 'callback function' when we press a key:
    switch(ev.keyCode){
        case 39:             //right arrow
            g_EyeX += 0.1;
            rayChange();
            draw(gl);
            break;
        case 37:             //left arrow
            g_EyeX -= 0.1;
            rayChange();
            draw(gl);
            break;
        case 38:             //up arrow
            g_EyeY += 0.1;
            rayChange();
            draw(gl);
            break;
        case 40:
            g_EyeY -= 0.1;   //down arrow
            rayChange();
            draw(gl);
            break;
        case 187:
            g_EyeZ += 0.1;   //+
            rayChange();
            draw(gl);
            break;
        case 189:
            g_EyeZ -= 0.1;   //-
            rayChange();
            draw(gl);
            break;
        case 65:             //a
            var xx = g_EyeX-g_AtX;
            var xz = g_EyeZ-g_AtZ;
            var dt = Math.sqrt(xx*xx + xz*xz);
            theta += 2;
            g_AtX = g_EyeX - dt*Math.sin(theta/180*Math.PI);
            g_AtZ = g_EyeZ - dt*Math.cos(theta/180*Math.PI);
            rayChange();
            draw(gl);
            break;
        case 68:              //d
            var xx1 = g_EyeX-g_AtX;
            var xz1 = g_EyeZ-g_AtZ;
            var dt1 = Math.sqrt(xx1*xx1 + xz1*xz1);
            theta -= 2;
            g_AtX = g_EyeX - dt1*Math.sin(theta/180*Math.PI);
            g_AtZ = g_EyeZ - dt1*Math.cos(theta/180*Math.PI);
            rayChange();
            draw(gl);
            break;
        case 87:              //w
            g_AtY += 0.1;
            rayChange();
            draw(gl);
            break;
        case 83:              //s
            g_AtY -= 0.1;
            rayChange();
            draw(gl);
            break;
        case 112:              //f1
            alert("User Instructions:"+ '\n'+ "(also brief instructions show on the canvas bellow)"+ '\n\n'+
                  "1. Press 't' or 'T' to see the ray tracing result shows on the right." + '\n'+
                  "2. Press 'n' to change the scene." + '\n'+
                  "3. Press 'l' to turn on/off the head light; 'k' to turn on/off the light one the eye."+ '\n'+
                  "4. Press 'p' to show aliasing/anti-aliasing result."+ '\n'+
                  "5. Press 'o' to add/remove the materials on the objects."+ '\n'+
                  "6. Press 'm' to add/remove shadow effect."+ '\n'+
                  "7. Press the up, down, left, right arrow key and '+' '-' to move the eye point."+ '\n'+
                  "8. Press the 'w', 's', 'a', 'd' to move the lookat point."+ '\n'+
                  "9. Click the lamp position buttons to change the head lamp 3D position."+ '\n'+
                  "10. Click the reflect++/reflect-- buttons to change the recursion depth for rays."+ '\n'+
                  "11. Move the mouse to transform the whole objects.");
            break;

    }
    
}

function lampXadd(){
    lampX += 0.1;
    draw(gl);
}
function lampXminus(){
    lampX -= 0.1;
    draw(gl);
}
function lampYadd(){
    lampY += 0.1;
    draw(gl);
}
function lampYminus(){
    lampY -= 0.1;
    draw(gl);
}
function lampZadd(){
    lampZ += 0.1;
    draw(gl);
}
function lampZminus(){
    lampZ -= 0.1;
    draw(gl);
}
function refDepthAdd(){
    myScene.depthMax += 1;
    document.getElementById('Result').innerHTML =
    ' Reflection Depth: ' + myScene.depthMax;
}
function refDepthMinus(){
    if (myScene.depthMax > 0){
        myScene.depthMax -= 1;
        document.getElementById('Result').innerHTML =
        ' Reflection Depth:' + myScene.depthMax;
    }
    else{
        myScene.depthMax = 0;
        document.getElementById('Result').innerHTML =
        ' Reflection Depth: ' + myScene.depthMax;
    }
    
}

function browserResize() {
    //==============================================================================
    // Called when user re-sizes their browser window , because our HTML file
    // contains:  <body onload="main()" onresize="browserResize()">
    
    /* SOLUTION to a pesky problem:
     The main() function retrieves our WebGL drawing context as the variable 'gl', then shares it as an argument to other functions.
     That's not enough!
     How can we access the 'gl' canvas within functions that main() will NEVER call, such as the mouse and keyboard-handling functions, or winResize()? Easy! make our own local references to the current canvas and WebGL drawing
     context, like this: */
    
    var myCanvas = document.getElementById('webgl');	// get current canvas
    var myGL = getWebGLContext(myCanvas);							// and context:
    //Report our current browser-window contents:
    
    //console.log('myCanvas width,height=', myCanvas.width, myCanvas.height);
    //    console.log('Browser window: innerWidth,innerHeight=',
    //																innerWidth, innerHeight);	// http://www.w3schools.com/jsref/obj_window.asp
    //
    //Make a square canvas/CVV fill the SMALLER of the width/2 or height:
    if(innerWidth > 2*innerHeight) {  // fit to brower-window height
        myCanvas.width = 2*innerHeight-20;
        myCanvas.height = innerHeight-20;
    }
    else {	// fit canvas to browser-window width
        myCanvas.width = innerWidth-20;
        myCanvas.height = 0.5*innerWidth-20;
    }
    //console.log('NEW myCanvas width,height=', myCanvas.width, myCanvas.height);
}

function makeGroundGrid() {
    //==============================================================================
    // Create a list of vertices that create a large grid of lines in the x,y plane
    // centered at x=y=z=0.  Draw this shape using the GL_LINES primitive.
    
    //    var xcount = 10;			// # of lines to draw in x,y to make the grid.
    //    var ycount = 10;
    //    var xymax	= 2.3;			// grid size; extends to cover +/-xymax in x and y.
    var xColr = new Float32Array([0.18, 0.85, 0.95]);	// bright yellow
    var yColr = new Float32Array([0.18, 0.85, 0.95]);	// bright green.
    
    // Create an (global) array to hold this ground-plane's vertices:
    gndVerts = new Float32Array(floatsPerVertex*2*(xcount+ycount));
    // draw a grid made of xcount+ycount lines; 2 vertices per line.
    
    var xgap = xymax/(xcount-1);		// HALF-spacing between lines in x,y;
    var ygap = xymax/(ycount-1);		// (why half? because v==(0line number/2))
    
    // First, step thru x values as we make vertical lines of constant-x:
    for(v=0, j=0; v<2*xcount; v++, j+= floatsPerVertex) {
        if(v%2==0) {	// put even-numbered vertices at (xnow, -xymax, 0)
            gndVerts[j  ] = -xymax + (v  )*xgap;	// x
            gndVerts[j+1] = -xymax;								// y
            gndVerts[j+2] = 0.0;									// z
            gndVerts[j+6] = -xymax + (v  )*xgap;
            gndVerts[j+7] = -xymax;
        }
        else {				// put odd-numbered vertices at (xnow, +xymax, 0).
            gndVerts[j  ] = -xymax + (v-1)*xgap;	// x
            gndVerts[j+1] = xymax;								// y
            gndVerts[j+2] = 0.0;									// z
            gndVerts[j+6] = -xymax + (v-1)*xgap;
            gndVerts[j+7] = xymax;
        }
        gndVerts[j+3] = xColr[0];			// red
        gndVerts[j+4] = xColr[1];			// grn
        gndVerts[j+5] = xColr[2];			// blu
        gndVerts[j+8] = 0.0;
        gndVerts[j+9] = 0.0;
        gndVerts[j+10] = 1.0;
    }
    // Second, step thru y values as wqe make horizontal lines of constant-y:
    // (don't re-initialize j--we're adding more vertices to the array)
    for(v=0; v<2*ycount; v++, j+= floatsPerVertex) {
        if(v%2==0) {		// put even-numbered vertices at (-xymax, ynow, 0)
            gndVerts[j  ] = -xymax;								// x
            gndVerts[j+1] = -xymax + (v  )*ygap;	// y
            gndVerts[j+2] = 0.0;									// z
            gndVerts[j+6] = -xymax;
            gndVerts[j+7] = -xymax + (v  )*ygap;;
        }
        else {					// put odd-numbered vertices at (+xymax, ynow, 0).
            gndVerts[j  ] = xymax;								// x
            gndVerts[j+1] = -xymax + (v-1)*ygap;	// y
            gndVerts[j+2] = 0.0;									// z
            gndVerts[j+6] = xymax;
            gndVerts[j+7] = -xymax + (v-1)*ygap;
        }
        gndVerts[j+3] = yColr[0];			// red
        gndVerts[j+4] = yColr[1];			// grn
        gndVerts[j+5] = yColr[2];			// blu
        gndVerts[j+8] = 0.0;
        gndVerts[j+9] = 0.0;
        gndVerts[j+10] = 1.0;
    }
}

function makeSphere(sColr) {
    //==============================================================================
    // Make a sphere from one OpenGL TRIANGLE_STRIP primitive.   Make ring-like
    // equal-lattitude 'slices' of the sphere (bounded by planes of constant z),
    // and connect them as a 'stepped spiral' design (see makeCylinder) to build the
    // sphere from one triangle strip.
    var slices = 30;		// # of slices of the sphere along the z axis. >=3 req'd
    // (choose odd # or prime# to avoid accidental symmetry)
    var sliceVerts	= 50;	// # of vertices around the top edge of the slice
    // (same number of vertices on bottom of slice, too)
    var topColr = new Float32Array([sColr[0], sColr[1], sColr[2]]);	// North Pole: light gray
    var equColr = new Float32Array([sColr[0], sColr[1], sColr[2]]);	// Equator:    bright green
    var botColr = new Float32Array([sColr[0], sColr[1], sColr[2]]);	// South Pole: brightest gray.
    var sliceAngle = Math.PI/slices;	// lattitude angle spanned by one slice.
    
    // Create a (global) array to hold this sphere's vertices:
    sphVerts = new Float32Array(  ((slices * 2* sliceVerts) -2) * floatsPerVertex);
    // # of vertices * # of elements needed to store them.
    // each slice requires 2*sliceVerts vertices except 1st and
    // last ones, which require only 2*sliceVerts-1.
    
    // Create dome-shaped top slice of sphere at z=+1
    // s counts slices; v counts vertices;
    // j counts array elements (vertices * elements per vertex)
    var cos0 = 0.0;					// sines,cosines of slice's top, bottom edge.
    var sin0 = 0.0;
    var cos1 = 0.0;
    var sin1 = 0.0;
    var j = 0;							// initialize our array index
    var isLast = 0;
    var isFirst = 1;
    for(s=0; s<slices; s++) {	// for each slice of the sphere,
        // find sines & cosines for top and bottom of this slice
        if(s==0) {
            isFirst = 1;	// skip 1st vertex of 1st slice.
            cos0 = 1.0; 	// initialize: start at north pole.
            sin0 = 0.0;
        }
        else {					// otherwise, new top edge == old bottom edge
            isFirst = 0;
            cos0 = cos1;
            sin0 = sin1;
        }								// & compute sine,cosine for new bottom edge.
        cos1 = Math.cos((s+1)*sliceAngle);
        sin1 = Math.sin((s+1)*sliceAngle);
        // go around the entire slice, generating TRIANGLE_STRIP verts
        // (Note we don't initialize j; grows with each new attrib,vertex, and slice)
        if(s==slices-1) isLast=1;	// skip last vertex of last slice.
        for(v=isFirst; v< 2*sliceVerts-isLast; v++, j+=floatsPerVertex) {
            if(v%2==0)
            {				// put even# vertices at the the slice's top edge
                // (why PI and not 2*PI? because 0 <= v < 2*sliceVerts
                // and thus we can simplify cos(2*PI(v/2*sliceVerts))
                sphVerts[j  ] = sin0 * Math.cos(Math.PI*(v)/sliceVerts);
                sphVerts[j+1] = sin0 * Math.sin(Math.PI*(v)/sliceVerts);
                sphVerts[j+2] = cos0;
                sphVerts[j+6] = sin0 * Math.cos(Math.PI*(v)/sliceVerts);
                sphVerts[j+7] = sin0 * Math.sin(Math.PI*(v)/sliceVerts);
                sphVerts[j+8] = sphVerts[j  ];
                sphVerts[j+9] = sphVerts[j+1];
                sphVerts[j+10] = sphVerts[j+2];
                
            }
            else { 	// put odd# vertices around the slice's lower edge;
                // x,y,z,w == cos(theta),sin(theta), 1.0, 1.0
                // 					theta = 2*PI*((v-1)/2)/capVerts = PI*(v-1)/capVerts
                sphVerts[j  ] = sin1 * Math.cos(Math.PI*(v-1)/sliceVerts);		// x
                sphVerts[j+1] = sin1 * Math.sin(Math.PI*(v-1)/sliceVerts);		// y
                sphVerts[j+2] = cos1;																				// z
                // w.
                sphVerts[j+6] = sin1 * Math.cos(Math.PI*(v-1)/sliceVerts);
                sphVerts[j+7] = sin1 * Math.sin(Math.PI*(v-1)/sliceVerts);
                sphVerts[j+8] = sphVerts[j  ];
                sphVerts[j+9] = sphVerts[j+1];
                sphVerts[j+10] = sphVerts[j+2];
            }
            if(s==0) {	// finally, set some interesting colors for vertices:
                sphVerts[j+3]=topColr[0];
                sphVerts[j+4]=topColr[1];
                sphVerts[j+5]=topColr[2];
            }
            else if(s==slices-1) {
                sphVerts[j+3]=botColr[0];
                sphVerts[j+4]=botColr[1];
                sphVerts[j+5]=botColr[2];
            }
            else {
                sphVerts[j+3]=botColr[0];// equColr[0];
                sphVerts[j+4]=botColr[1];// equColr[1];
                sphVerts[j+5]=botColr[2];// equColr[2];
            }
            
        }
    }
}

function makeCube(){
    cube = new Float32Array([
      -0.5,  0.5,  0.5,  0.6, 0.97, 0.55,  -0.5, 0.5,  0.0, 0.0, 1.0,       //front
      -0.5, -0.5,  0.5,  0.6, 0.97, 0.55,  -0.5,-0.5,  0.0, 0.0, 1.0,
       0.5,  0.5,  0.5,  0.6, 0.97, 0.55,   0.5, 0.5,  0.0, 0.0, 1.0,
       
      -0.5, -0.5,  0.5,  0.6, 0.97, 0.55,  -0.5,-0.5,  0.0, 0.0, 1.0,
       0.5,  0.5,  0.5,  0.6, 0.97, 0.55,   0.5, 0.5,  0.0, 0.0, 1.0,
       0.5, -0.5,  0.5,  0.6, 0.97, 0.55,   0.5,-0.5,  0.0, 0.0, 1.0,
         
       0.5,  0.5,  0.5,  0.6, 0.97, 0.55,   0.5, 0.5,  1.0, 0.0, 0.0,       //right
       0.5, -0.5,  0.5,  0.6, 0.97, 0.55,   0.5,-0.5,  1.0, 0.0, 0.0,
       0.5,  0.5, -0.5,  0.6, 0.97, 0.55,   0.5, 0.5,  1.0, 0.0, 0.0,
                             
       0.5, -0.5,  0.5,  0.6, 0.97, 0.55,   0.5,-0.5,  1.0, 0.0, 0.0,
       0.5,  0.5, -0.5,  0.6, 0.97, 0.55,   0.5, 0.5,  1.0, 0.0, 0.0,
       0.5, -0.5, -0.5,  0.6, 0.97, 0.55,   0.5,-0.5,  1.0, 0.0, 0.0,
       
       0.5,  0.5, -0.5,  0.6, 0.97, 0.55,   0.5, 0.5,  0.0, 0.0,-1.0,       //back
       0.5, -0.5, -0.5,  0.6, 0.97, 0.55,   0.5,-0.5,  0.0, 0.0,-1.0,
      -0.5,  0.5, -0.5,  0.6, 0.97, 0.55,  -0.5, 0.5,  0.0, 0.0,-1.0,
                             
       0.5, -0.5, -0.5,  0.6, 0.97, 0.55,   0.5,-0.5,  0.0, 0.0,-1.0,
      -0.5,  0.5, -0.5,  0.6, 0.97, 0.55,  -0.5, 0.5,  0.0, 0.0,-1.0,
      -0.5, -0.5, -0.5,  0.6, 0.97, 0.55,  -0.5,-0.5,  0.0, 0.0,-1.0,
                             
      -0.5,  0.5, -0.5,  0.6, 0.97, 0.55,  -0.5, 0.5, -1.0, 0.0, 0.0,       //left
      -0.5, -0.5, -0.5,  0.6, 0.97, 0.55,  -0.5,-0.5, -1.0, 0.0, 0.0,
      -0.5,  0.5,  0.5,  0.6, 0.97, 0.55,  -0.5, 0.5, -1.0, 0.0, 0.0,
                             
      -0.5, -0.5, -0.5,  0.6, 0.97, 0.55,  -0.5,-0.5, -1.0, 0.0, 0.0,
      -0.5,  0.5,  0.5,  0.6, 0.97, 0.55,  -0.5, 0.5, -1.0, 0.0, 0.0,
      -0.5, -0.5,  0.5,  0.6, 0.97, 0.55,  -0.5,-0.5, -1.0, 0.0, 0.0,
                             
      -0.5,  0.5, -0.5,  0.6, 0.97, 0.55,  -0.5, 0.5,  0.0, 1.0, 0.0,       //up
      -0.5,  0.5,  0.5,  0.6, 0.97, 0.55,  -0.5, 0.5,  0.0, 1.0, 0.0,
       0.5,  0.5, -0.5,  0.6, 0.97, 0.55,   0.5, 0.5,  0.0, 1.0, 0.0,
                             
      -0.5,  0.5,  0.5,  0.6, 0.97, 0.55,  -0.5, 0.5,  0.0, 1.0, 0.0,
       0.5,  0.5, -0.5,  0.6, 0.97, 0.55,   0.5, 0.5,  0.0, 1.0, 0.0,
       0.5,  0.5,  0.5,  0.6, 0.97, 0.55,   0.5, 0.5,  0.0, 1.0, 0.0,
                             
      -0.5, -0.5, -0.5,  0.6, 0.97, 0.55,  -0.5,-0.5,  0.0,-1.0, 0.0,       //bottom
      -0.5, -0.5,  0.5,  0.6, 0.97, 0.55,  -0.5,-0.5,  0.0,-1.0, 0.0,
       0.5, -0.5, -0.5,  0.6, 0.97, 0.55,   0.5,-0.5,  0.0,-1.0, 0.0,
                             
      -0.5, -0.5,  0.5,  0.6, 0.97, 0.55,  -0.5,-0.5,  0.0,-1.0, 0.0,
       0.5, -0.5, -0.5,  0.6, 0.97, 0.55,   0.5,-0.5,  0.0,-1.0, 0.0,
       0.5, -0.5,  0.5,  0.6, 0.97, 0.55,   0.5,-0.5,  0.0,-1.0, 0.0,
                             
    ]);
}

function makeCylinder(tR){
    //==============================================================================
    // Make a cylinder shape from one TRIANGLE_STRIP drawing primitive, using the
    // 'stepped spiral' design described in notes.
    // Cylinder center at origin, encircles z axis, radius 1, top/bottom at z= +/-1.
    //
    var ctrColr = new Float32Array([87/255, 96/255, 105/255]);	// dark gray
    var topColr = new Float32Array([249/255, 205/255, 173/255]);	// light green
    var botColr = new Float32Array([249/255, 205/255, 173/255]);	// light blue
    var capVerts = 40;	// # of vertices around the topmost 'cap' of the shape
    var botRadius = 1.0;		// radius of bottom of cylinder (top always 1.0)
    var topRadius = tR;
    
    // Create a (global) array to hold this cylinder's vertices;
    cylVerts = new Float32Array(  (capVerts*6+6)  * floatsPerVertex);
    // # of vertices * # of elements needed to store them.
    
    // Create circle-shaped top cap of cylinder at z=+1.0, radius 1.0
    // v counts vertices: j counts array elements (vertices * elements per vertex)
    for(v=1,j=0; v<2*capVerts+3; v++,j+=floatsPerVertex) {
        // skip the first vertex--not needed.
        if(v%2==0)
        {				// put even# vertices at center of cylinder's top cap:
            cylVerts[j  ] = 0.0; 			// x,y,z,w == 0,0,1,1
            cylVerts[j+1] = 0.0;
            cylVerts[j+2] = 1.0;
            cylVerts[j+3]= ctrColr[0];
            cylVerts[j+4]= ctrColr[1];
            cylVerts[j+5]= ctrColr[2];
            cylVerts[j+6]= 0.0;
            cylVerts[j+7]= 0.0;
            cylVerts[j+8]= 0.0;
            cylVerts[j+9]= 0.0;
            cylVerts[j+10]= 1.0;
        }
        else { 	// put odd# vertices around the top cap's outer edge;
            // x,y,z,w == cos(theta),sin(theta), 1.0, 1.0
            // 					theta = 2*PI*((v-1)/2)/capVerts = PI*(v-1)/capVerts
            cylVerts[j  ] = tR * Math.cos(Math.PI*(v-1)/capVerts);			// x
            cylVerts[j+1] = tR * Math.sin(Math.PI*(v-1)/capVerts);			// y
            //	(Why not 2*PI? because 0 < =v < 2*capVerts, so we
            //	 can simplify cos(2*PI * (v-1)/(2*capVerts))
            cylVerts[j+2] = 1.0;	// z
            // r,g,b = topColr[]
            cylVerts[j+3]= topColr[0];
            cylVerts[j+4]= topColr[1];
            cylVerts[j+5]= topColr[2];
            cylVerts[j+6]= Math.cos(Math.PI*(v-1)/capVerts);
            cylVerts[j+7]= Math.sin(Math.PI*(v-1)/capVerts);
            cylVerts[j+8]= 0.0;
            cylVerts[j+9]= 0.0;
            cylVerts[j+10]= 1.0;
        }
    }
    // Create the cylinder side walls, made of 2*capVerts vertices.
    // v counts vertices within the wall; j continues to count array elements
    for(v=0; v< 2*capVerts+2; v++, j+=floatsPerVertex) {
        if(v%2==0)	// position all even# vertices along top cap:
        {
            cylVerts[j  ] = tR * Math.cos(Math.PI*(v)/capVerts);		// x
            cylVerts[j+1] = tR * Math.sin(Math.PI*(v)/capVerts);		// y
            cylVerts[j+2] = 1.0;	// z
            // r,g,b = topColr[]
            cylVerts[j+3]= topColr[0];
            cylVerts[j+4]= topColr[1];
            cylVerts[j+5]= topColr[2];
            cylVerts[j+6]= Math.cos(Math.PI*(v)/capVerts);
            cylVerts[j+7]= Math.sin(Math.PI*(v)/capVerts);
            cylVerts[j+8]= Math.cos(Math.PI*(v)/capVerts);
            cylVerts[j+9]= Math.sin(Math.PI*(v)/capVerts);
            cylVerts[j+10]= -(tR-1)*(1+(tR-1));
        }
        else		// position all odd# vertices along the bottom cap:
        {
            cylVerts[j  ] = botRadius * Math.cos(Math.PI*(v-1)/capVerts);		// x
            cylVerts[j+1] = botRadius * Math.sin(Math.PI*(v-1)/capVerts);		// y
            cylVerts[j+2] = 0.0;	// z
            // r,g,b = topColr[]
            cylVerts[j+3]= botColr[0];
            cylVerts[j+4]= botColr[1];
            cylVerts[j+5]= botColr[2];
            cylVerts[j+6]= botRadius * Math.cos(Math.PI*(v-1)/capVerts);
            cylVerts[j+7]= botRadius * Math.sin(Math.PI*(v-1)/capVerts);
            cylVerts[j+8]= botRadius * Math.cos(Math.PI*(v-1)/capVerts);
            cylVerts[j+9]= botRadius * Math.sin(Math.PI*(v-1)/capVerts);
            cylVerts[j+10]= -(tR-1)*(1-(tR-1));
        }
    }
    // Create the cylinder bottom cap, made of 2*capVerts -1 vertices.
    // v counts the vertices in the cap; j continues to count array elements
    for(v=0; v < 2*capVerts+1; v++, j+= floatsPerVertex) {
        if(v%2==0) {	// position even #'d vertices around bot cap's outer edge
            cylVerts[j  ] = botRadius * Math.cos(Math.PI*(v)/capVerts);		// x
            cylVerts[j+1] = botRadius * Math.sin(Math.PI*(v)/capVerts);		// y
            cylVerts[j+2] = 0.0;	// z
            // r,g,b = topColr[]
            cylVerts[j+3]=botColr[0];
            cylVerts[j+4]=botColr[1];
            cylVerts[j+5]=botColr[2];
            cylVerts[j+6]= botRadius * Math.cos(Math.PI*(v)/capVerts);
            cylVerts[j+7]= botRadius * Math.sin(Math.PI*(v)/capVerts);
            cylVerts[j+8]= 0.0;
            cylVerts[j+9]= 0.0;
            cylVerts[j+10]= -1.0;
        }
        else {				// position odd#'d vertices at center of the bottom cap:
            cylVerts[j  ] = 0.0; 			// x,y,z,w == 0,0,-1,1
            cylVerts[j+1] = 0.0;	
            cylVerts[j+2] = 0.0;
            cylVerts[j+3]=botColr[0];
            cylVerts[j+4]=botColr[1];
            cylVerts[j+5]=botColr[2];
            cylVerts[j+6]= 0.0;
            cylVerts[j+7]= 0.0;
            cylVerts[j+8]= 0.0;
            cylVerts[j+9]= 0.0;
            cylVerts[j+10]=-1.0;
        }
    }
}

function dragQuat(xdrag, ydrag) {
    //==============================================================================
    // Called when user drags mouse by 'xdrag,ydrag' as measured in CVV coords.
    // We find a rotation axis perpendicular to the drag direction, and convert the
    // drag distance to an angular rotation amount, and use both to set the value of
    // the quaternion qNew.  We then combine this new rotation with the current
    // rotation stored in quaternion 'qTot' by quaternion multiply.  Note the
    // 'draw()' function converts this current 'qTot' quaternion to a rotation
    // matrix for drawing.
    var res = 5;
    var qTmp = new Quaternion(0,0,0,1);
    
    var dist = Math.sqrt(xdrag*xdrag + ydrag*ydrag);
    // console.log('xdrag,ydrag=',xdrag.toFixed(5),ydrag.toFixed(5),'dist=',dist.toFixed(5));
    qNew.setFromAxisAngle(-ydrag + 0.0001, xdrag + 0.0001, 0.0, dist*150.0);
    // (why add tiny 0.0001? To ensure we never have a zero-length rotation axis)
    // why axis (x,y,z) = (-yMdrag,+xMdrag,0)?
    // -- to rotate around +x axis, drag mouse in -y direction.
    // -- to rotate around +y axis, drag mouse in +x direction.
    qTmp.multiply(qNew,qTot);			// apply new rotation to current rotation.
    qTot.copy(qTmp);
    //--------------------------
    // IMPORTANT! Why qNew*qTot instead of qTot*qNew? (Try it!)
    // ANSWER: Because 'duality' governs ALL transformations, not just matrices.
    // If we multiplied in (qTot*qNew) order, we would rotate the drawing axes
    // first by qTot, and then by qNew--we would apply mouse-dragging rotations
    // to already-rotated drawing axes.  Instead, we wish to apply the mouse-drag
    // rotations FIRST, before we apply rotations from all the previous dragging.
    //------------------------
    // IMPORTANT!  Both qTot and qNew are unit-length quaternions, but we store
    // them with finite precision. While the product of two (EXACTLY) unit-length
    // quaternions will always be another unit-length quaternion, the qTmp length
    // may drift away from 1.0 if we repeat this quaternion multiply many times.
    // A non-unit-length quaternion won't work with our quaternion-to-matrix fcn.
    // Matrix4.prototype.setFromQuat().
    //	qTmp.normalize();						// normalize to ensure we stay at length==1.0.
    
    // show the new quaternion qTot on our webpage in the <div> element 'QuatValue'
//    document.getElementById('QuatValue').innerHTML=
//    '\t X=' +qTot.x.toFixed(res)+
//    'i\t Y=' +qTot.y.toFixed(res)+
//    'j\t Z=' +qTot.z.toFixed(res)+
//    'k\t W=' +qTot.w.toFixed(res)+
//    '<br>length='+qTot.length().toFixed(res);
};

function testQuaternions() {
    //==============================================================================
    // Test our little "quaternion-mod.js" library with simple rotations for which
    // we know the answers; print results to make sure all functions work as
    // intended.
    // 1)  Test constructors and value-setting functions:
    
    var res = 5;
    var myQuat = new Quaternion(1,2,3,4);
    console.log('constructor: myQuat(x,y,z,w)=',
                myQuat.x, myQuat.y, myQuat.z, myQuat.w);
    myQuat.clear();
    console.log('myQuat.clear()=',
                myQuat.x.toFixed(res), myQuat.y.toFixed(res),
                myQuat.z.toFixed(res), myQuat.w.toFixed(res));
    myQuat.set(1,2, 3,4);
    console.log('myQuat.set(1,2,3,4)=',
                myQuat.x.toFixed(res), myQuat.y.toFixed(res),
                myQuat.z.toFixed(res), myQuat.w.toFixed(res));
    console.log('myQuat.length()=', myQuat.length().toFixed(res));
    myQuat.normalize();
    console.log('myQuat.normalize()=',
                myQuat.x.toFixed(res), myQuat.y.toFixed(res), myQuat.z.toFixed(res), myQuat.w.toFixed(res));
    // Simplest possible quaternions:
    myQuat.setFromAxisAngle(1,0,0,0);
    console.log('Set myQuat to 0-deg. rot. on x axis=',
                myQuat.x.toFixed(res), myQuat.y.toFixed(res), myQuat.z.toFixed(res), myQuat.w.toFixed(res));
    myQuat.setFromAxisAngle(0,1,0,0);
    console.log('set myQuat to 0-deg. rot. on y axis=',
                myQuat.x.toFixed(res), myQuat.y.toFixed(res), myQuat.z.toFixed(res), myQuat.w.toFixed(res));
    myQuat.setFromAxisAngle(0,0,1,0);
    console.log('set myQuat to 0-deg. rot. on z axis=',
                myQuat.x.toFixed(res), myQuat.y.toFixed(res), myQuat.z.toFixed(res), myQuat.w.toFixed(res), '\n');
    
    myQmat = new Matrix4();
    myQuat.setFromAxisAngle(1,0,0, 90.0);
    console.log('set myQuat to +90-deg rot. on x axis =',
                myQuat.x.toFixed(res), myQuat.y.toFixed(res), myQuat.z.toFixed(res), myQuat.w.toFixed(res));
    myQmat.setFromQuat(myQuat.x, myQuat.y, myQuat.z, myQuat.w);
    console.log('myQuat as matrix: (+y axis <== -z axis)(+z axis <== +y axis)');
    myQmat.printMe();
    
    myQuat.setFromAxisAngle(0,1,0, 90.0);
    console.log('set myQuat to +90-deg rot. on y axis =',
                myQuat.x.toFixed(res), myQuat.y.toFixed(res), myQuat.z.toFixed(res), myQuat.w.toFixed(res));
    myQmat.setFromQuat(myQuat.x, myQuat.y, myQuat.z, myQuat.w);
    console.log('myQuat as matrix: (+x axis <== +z axis)(+z axis <== -x axis)');
    myQmat.printMe();
    
    myQuat.setFromAxisAngle(0,0,1, 90.0);
    console.log('set myQuat to +90-deg rot. on z axis =',
                myQuat.x.toFixed(res), myQuat.y.toFixed(res), myQuat.z.toFixed(res), myQuat.w.toFixed(res));
    myQmat.setFromQuat(myQuat.x, myQuat.y, myQuat.z, myQuat.w);
    console.log('myQuat as matrix: (+x axis <== -y axis)(+y axis <== +x axis)');
    myQmat.printMe();
    
    // Test quaternion multiply:
    // (q1*q2) should rotate drawing axes by q1 and then by q2;  it does!
    var qx90 = new Quaternion;
    var qy90 = new Quaternion;
    qx90.setFromAxisAngle(1,0,0,90.0);			// +90 deg on x axis
    qy90.setFromAxisAngle(0,1,0,90.0);			// +90 deg on y axis.
    myQuat.multiply(qx90,qy90);
    console.log('set myQuat to (90deg x axis) * (90deg y axis) = ',
                myQuat.x.toFixed(res), myQuat.y.toFixed(res), myQuat.z.toFixed(res), myQuat.w.toFixed(res));
    myQmat.setFromQuat(myQuat.x, myQuat.y, myQuat.z, myQuat.w);
    console.log('myQuat as matrix: (+x <== +z)(+y <== +x )(+z <== +y');
    myQmat.printMe();
}
