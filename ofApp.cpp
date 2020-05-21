#include "ofApp.h"
float G = 6.67e-11;
int n_max = 1; float dt = 10000; float theta = 0.3; 
float m_n = 5000; float r1 = 1000; float M1 = 60;

class Axis {
public:
	float m_x11; float m_x12; float m_x21; float m_x22;
	float m_y11; float m_y12; float m_y21; float m_y22;

	void set(float x11, float x12, float x21, float x22,
		float y11, float y12, float y21, float y22) {
		m_x11 = x11; m_x12 = x12; m_x21 = x21; m_x22 = x22;
		m_y11 = y11; m_y12 = y12; m_y21 = y21; m_y22 = y22;
	}

	void show() {
		int n = 1;

		float dx = (m_x12 / m_x22+30) / n;
		int n1 = floor(m_x12 / dx);

		float dy = (m_y12 / m_y22+30) / n;
		int n2 = floor(m_y12 / dy);

		ofSetColor(220, 220, 220, 10);

		for (int i = 0; i < n1; i++) {
			float x1 = i * dx;
			float x2 = -(i + 1)*dx;
			ofDrawLine(x1, m_y11, x1, m_y12);
			ofDrawLine(x2, m_y11, x2, m_y12);
		}

		for (int i = 0; i < n2; i++) {
			float y1 = i * dx;
			float y2 = -(i + 1)*dx;
			ofDrawLine(m_x11, y1, m_x12, y1);
			ofDrawLine(m_x11, y2, m_x12, y2);
		}
	}

	void e(float x11, float x21, float r11) {
		float x12 = ofMap(x11, m_x21, m_x22, m_x11, m_x12);
		float x22 = ofMap(x21, m_y21, m_y22, m_y11, m_y12);
		float r12 = ofMap(r11, 0, m_x22, 0, m_x12);
		ofDrawCircle(x12, x22, r11);
	}

	void r(float x11, float y11, float w1, float h1) {
		ofSetRectMode(OF_RECTMODE_CENTER);
		float x12 = ofMap(x11, m_x21, m_x22, m_x11, m_x12);
		float y12 = ofMap(y11, m_y21, m_y22, m_y11, m_y12);
		float w2 = ofMap(w1, 0, m_x22, 0, 2*m_x12);
		float h2 = ofMap(h1, 0, m_y22, 0, 2*m_y12);
		ofDrawRectangle(x12, y12, w2, h2);
	}

	void l(float x11, float y11, float x21, float y21) {
		float x12 = ofMap(x11, m_x21, m_x22, m_x11, m_x12);
		float y12 = ofMap(y11, m_y21, m_y22, m_y11, m_y12);
		float x22 = ofMap(x21, m_x21, m_x22, m_x11, m_x12);
		float y22 = ofMap(y21, m_y21, m_y22, m_y11, m_y12);
		ofDrawLine(x12, y12, x22, y22);
	}
};

ofVec2f Grav(ofVec2f x1, ofVec2f x2, float m2) {
	ofVec2f x3; float d; float k;
	x3 = x2 - x1; d = pow(x3.length(), 2);
	k = pow(d+1,1.5);
	x3 *= (G*m2/k);
	return x3;
}

ofVec2f paralel(ofVec2f v1) {
	ofVec2f v2; float x11, x12;
	x11 = v1.x; x12 = v1.y; v2.set(-x12, x11);
	return v2;
}

ofVec2f direction(ofVec2f v1, float mag) {
	ofVec2f v2; float x11, x12, x21, x22, s, c;
	v1.normalize();
	x11 = v1.x; x12 = v1.y; s = 1;

	if (x11 != 0) {
		if (x11 >= 0) {
			s = -1;
		}
		c = x12 / x11;
		x21 = -s * pow((pow(mag, 2)) / ((pow(c, 2) + 1)), .5);
		x22 = c * x21;
		v2.set(x21, x22);
		return v2;
	}
	if (x12 != 0) {
		if (x12 >= 0) {
			s = -1;
		}
		c = x11 / x12;
		x22 = -s * pow((pow(mag, 2)) / ((pow(c, 2) + 1)), .5);
		x21 = c * x22;
		v2.set(x21, x22);
		return v2;
	}
	if (x12 == 0 && x11 == 0) {
		v2.set(0, 0);
		return v2;
	}
}

vector<ofVec2f> distribution(int n, float r) {
	vector<ofVec2f> v1;
	for (int i = 0; i < n; i++) {
		float x1, x2; ofVec2f v2; 
		x1 = ofRandom(-r, r); x2 = ofRandom(-r, r);
		if (x2 <= 0 && x2 >= -(pow((pow(r, 2) - pow(x1, 2)), .5)))
		{
			v2.set(x1, x2);
			v1.push_back(v2);
		}
		if (x2 >= 0 && x2 <= pow((pow(r, 2) - pow(x1, 2)), .5)) {
			v2.set(x1, x2);
			v1.push_back(v2);
		}
	}
	return v1;
}

class Particle {
public:
	ofVec2f m_x, m_v, m_a, m_a2;
	float mass;

	void set(ofVec2f x, ofVec2f v, ofVec2f a, float m) {
		m_x = x; m_v = v; m_a = a; mass = m; m_a2.set(0, 0);
	}

	void show(Axis& a) {
		a.e(m_x.x, m_x.y, 1);
	}
};

void up1(Particle * p1, ofVec2f * x21, float * m2) {
	ofVec2f a11;  
	a11 = Grav(p1->m_x, *x21, *m2);
	p1->m_a += a11;
}

void up2(Particle * p1, ofVec2f * x22, float * m2) {
	ofVec2f * x12; ofVec2f * v11; ofVec2f v12; ofVec2f * a11;
	x12 = &p1->m_x; v11 = &p1->m_v; a11 = &p1->m_a;
	p1->m_a2 += Grav(*x12, *x22, *m2);
}

class Particles {
public:
	vector<Particle> P;

	void set(int n) {
		int m_n = n;
		for (int i = 0; i < m_n; i++)
		{
			Particle p;
			P.push_back(p);
		}
		vector<ofVec2f> pos; pos = distribution(m_n, r1);
		for (int i = 0; i < pos.size(); i++) {
			ofVec2f x, v1, v2, a; Particle p; float m, vo, mo;
			x = pos[i];
			//mo = exp(ofMap(x.length(), 0, r1, 50, 0));
			mo = ofRandom(1, 1000);
			//vo = pow((G*exp(M1)/x.length()), .5);
			//v1 = paralel(x);
			//v2 = direction(v1, vo);
			v2.set(0, 0);
			a.set(0, 0); m = mo;
			p.set(x, v2, a, m);
			P.push_back(p);
		}

		/*P[0].m_x.set(-1.0e11, 0); P[1].m_x.set(3.0e11, 0); P[2].m_x.set(5.0e11, 0);
		P[0].m_v.set(0,-9132.360046296); P[1].m_v.set(0, 27397.080138889); P[2].m_v.set(0, 47106.48148);
		P[0].mass = 6.0e30; P[1].mass = 2.0e30; P[2].mass = 1.0e27;*/
	}
	void update() {
		for (int i = 0; i < P.size(); i++) {
			ofVec2f *x11; ofVec2f x12; ofVec2f *v11; ofVec2f a11; float *m1;
			x11 = &P[i].m_x; v11 = &P[i].m_v; m1 = &P[i].mass;
			a11.set(0, 0);
			for (int j = 0; j < P.size(); j++) {
				ofVec2f* x21; float* m2;
				x21 = &P[j].m_x, m2 = &P[j].mass;
				if (i == j) {
					continue;
				}
				else {
					a11 += Grav(*x11, *x21, *m2);
				}
			}
			P[i].m_a = a11;
		}
		for (int i = 0; i < P.size(); i++) {
			ofVec2f * x11; ofVec2f x12; ofVec2f * v11; ofVec2f * a11; float * m1;
			x11 = &P[i].m_x; v11 = &P[i].m_v; m1 = &P[i].mass;
			a11 = &P[i].m_a;
			x12 = *x11 + *v11*dt + *a11*(.5*dt*dt);
			P[i].m_x = x12;
		}
		for (int i = 0; i < P.size(); i++) {
			ofVec2f * x12; ofVec2f * v11; ofVec2f v12; ofVec2f * a11; ofVec2f a12; float * m1;
			x12 = &P[i].m_x, v11 = &P[i].m_v, a11 = &P[i].m_a, m1 = &P[i].mass;
			a12.set(0.0);
			for (int j = 0; j < P.size(); j++) {
				ofVec2f * x22; float * m2;
				x22 = &P[j].m_x, m2 = &P[j].mass;
				if (i == j) {
					continue;
				}
				else {
					a12 += Grav(*x12, *x22, *m2);
				}
			}
			v12 = *v11 + (*a11 + a12)*(.5*dt);
			P[i].m_v = v12;
		}
	}
	void show(Axis& a) {
		for (int i = 0; i < P.size(); i++) {
			P[i].show(a);
		}
	}
};

class Rect {
public:
	float x, y, w, h;
	void set(float x1, float y1, float w1, float h1) {
		x = x1; y = y1;
		w = w1; h = h1;
	}
	bool contains(ofVec2f point) {
		bool b;
		b = false;
		if (x - w < point.x && point.x <= x + w
			&& y - h < point.y && point.y <= y + h ) {
			b = true;
		}
		return b;
	}
};

class Node {
public:
	Rect rect;
	Node * parent;
	vector<Node> C;
	vector<Particle *> P;
	ofVec2f CM; float mass;
	bool t1 = false; bool t2 = false; bool t3 = false; bool divided = false;

	void set(Rect r) {
		rect = r;
	}
	void show(Axis& a) {
		if (divided) {
			for (int i = 0; i < C.size(); i++) {
				C[i].show(a);
			}
		}
		a.r(rect.x, rect.y, rect.w, rect.h);
	}
	void type() {
		if (divided) {
			for (int i = 0; i < C.size(); i++) {
				C[i].type();
			}
		}
		if (P.size() == n_max) {
			t1 = true;
		}
		if (divided) {
			t2 = true;
		}
		if (C.size() == 0 && P.size() == 0) {
			t3 = true;
		}
	}
	void center() {
		if (divided = true) {
			for (int i = 0; i < C.size(); i++) {
				C[i].center();
			}
		}
		if (t1) {
			CM = P[0]->m_x;
			mass = P[0]->mass;
		}
		if (t2) {
			CM.set(0, 0);
			mass = 0;
			for (int i = 0; i < C.size(); i++) {
				CM += C[i].CM*C[i].mass;
				mass += C[i].mass;
			}
			CM *= (1/mass);
		}
		if (t3) {
			CM.set(0, 0);
			mass = 0;
		}
	}
	void subdivide() {
		Node nw, ne, sw, se;
		Rect nw_r, ne_r, sw_r, se_r;
		float * x; float * y; float * w; float * h;
		x = &rect.x; y = &rect.y; w = &rect.w; h = &rect.h; 
		nw_r.set(*x - *w / 2, *y + *h / 2, *w / 2, *h / 2); nw.set(nw_r);
		ne_r.set(*x + *w / 2, *y + *h / 2, *w / 2, *h / 2); ne.set(ne_r);
		sw_r.set(*x - *w / 2, *y - *h / 2, *w / 2, *h / 2); sw.set(sw_r);
		se_r.set(*x + *w / 2, *y - *h / 2, *w / 2, *h / 2); se.set(se_r);
		C.push_back(nw);
		C.push_back(ne);
		C.push_back(sw);
		C.push_back(se);
		for (int i = 0; i <= n_max; i++) {
			for (int j = 0; j < C.size(); j++) {
				C[j].insert(P[n_max-i]);
				C[j].parent = this;
			}
			P.pop_back();
		}
		divided = true;
	}
	void insert(Particle * planet) {
		if (rect.contains(planet->m_x) == false) {
			return;
		}
		if (divided) {
			for (int i = 0; i < C.size(); i++) {
				C[i].insert(planet);
			}
		}
		if (P.size() <= n_max && divided == false) {
			P.push_back(planet);
		}
		if (P.size() > n_max && divided == false) {
			subdivide();
		}
	}
	void interaction(Axis& a, Node * root) {
		if (divided) {
			for (int i = 0; i < C.size(); i++) {
				C[i].interaction(a, root);
			}
		}
		if (!t2) {
			return;
		}
		for (int i = 0; i < C.size(); i++) {
			if (!C[i].t1) {
				continue;
			}
			else {
				//cout << C[i].P[0]->m_x << "\n";
				Particle * p; p = C[i].P[0];
				ofVec2f * x11; ofVec2f * v11; ofVec2f a11; ofVec2f x12;
				x11 = &p->m_x; v11 = &p->m_v; a11 = p->m_a;
				x12 = *x11 + *v11*dt + a11 * (.5*dt*dt);
				p->m_x = x12;
			}
		}
	}
	void interaction2(Axis& a, Node * root) {
		if (divided) {
			for (int i = 0; i < C.size(); i++) {
				C[i].interaction2(a, root);
			}
		}
		if (!t2) {
			return;
		}
		for (int i = 0; i < C.size(); i++) {
			if (!C[i].t1) {
				continue;
			}
			else {
				reaches2(&C[i], root, a);
				Particle * p; p = C[i].P[0];
				ofVec2f * v11; ofVec2f * a11; ofVec2f * a12; ofVec2f v12;
				v11 = &p->m_v; a11 = &p->m_a; a12 = &p->m_a2; 
				v12 = *v11 + (*a11 + *a12)*(.5*dt);
				p->m_a = *a12;
				p->m_v = v12;
				p->m_a2.set(0, 0);
			}
		}
	}
	void reaches2(Node * n1, Node * n2, Axis& a) {
		float d, r; float * s; bool reach = false;
		d = (n2->CM - n1->CM).length();
		s = &n2->rect.w;
		r = *s / d;
		if (r < theta|| n2->t1) {
			up2(n1->P[0], &(*n2).CM, &(*n2).mass);
		}
		else {
			for (int i = 0; i < C.size(); i++) {
				if (n1 == &n2->C[i] || n2->C[i].t3) {
					continue;
				}
				else {
					reaches2(n1, &n2->C[i], a);
				}
			}
		}
	}
	void reaches(Node * n1, Node * n2, Axis& a) {
		float d, r; float * s; bool reach = false;
		d = (n2->CM - n1->CM).length();
		s = &n2->rect.w;
		r = *s / d;
		if (r < theta || n2->t1) {
			up1(n1->P[0], &(*n2).CM, &(*n2).mass);
		}
		else {
			for (int i = 0; i < C.size(); i++) {
				if (n1 == &n2->C[i] || n2->C[i].t3) {
					continue;
				}
				else {
					reaches(n1, &n2->C[i], a);
				}
			}
		}
	}
};

class quadTree {
public:
	Node root;

	void set(Rect bound) {
		root.set(bound);
	}
};

Axis a; Particles Planets; 

void resetSketch(float x11_i, float x12_i, float x21_i, float x22_i,
	float y11_i, float y12_i, float y21_i, float y22_i) {
	float x11 = x11_i;
	float x12 = x12_i;
	float x21 = x21_i;
	float x22 = x22_i;
	float y11 = y11_i;
	float y12 = y12_i;
	float y21 = y21_i;
	float y22 = y22_i;
	a.set(x11, x12, x21, x22, y11, y12, y21, y22);
}

//--------------------------------------------------------------
void ofApp::setup() {
	ofBackground(10);
	ofSetFrameRate(60);
	float x11 = -600.0;
	float x12 = 600.0;
	float x21 = -r1*2;
	float x22 = r1*2;
	float y11 = -300.0;
	float y12 = 300.0;
	float y21 = -(x22/x12)*y12;
	float y22 = (x22/x12)*y12;
	Planets.set(m_n);
	resetSketch(x11, x12, x21, x22, y11, y12, y21, y22);
}

//--------------------------------------------------------------
void ofApp::update() {

}
//--------------------------------------------------------------
void ofApp::draw() {
	ofTranslate(ofGetWidth() / 2, ofGetHeight() / 2); ofScale(1, -1, 1);
	ofNoFill();
	a.show();
	ofSetColor(200,30);
	quadTree Tree;
	Rect r; r.set(0, 0, r1, r1); Tree.set(r);
	for (int i = 0; i < Planets.P.size(); i++) {
		Tree.root.insert(&Planets.P[i]);
	}
	Tree.root.type();
	Tree.root.center();
	Tree.root.interaction(a, &Tree.root);
	Tree.root.interaction2(a, &Tree.root);
	Tree.root.show(a);
	ofSetColor(255);
	ofFill();
	Planets.show(a);
	cout << ofGetFrameRate() << "\n";
}


//--------------------------------------------------------------
void ofApp::keyPressed(int key) {
	float x; float y;
	if (key == OF_KEY_UP) {
		x = a.m_x22 + r1/10; y = (x / a.m_x12)*a.m_y12;
		resetSketch(a.m_x11, a.m_x12, -x, x, a.m_y11, a.m_y12, -y, y);
	}
	else if (key == OF_KEY_DOWN) {
		x = a.m_x22 - r1/10; y = (x / a.m_x12)*a.m_y12;
		resetSketch(a.m_x11, a.m_x12, -x, x, a.m_y11, a.m_y12, -y, y);
	}
	else if (key == OF_KEY_LEFT) {
		dt = 1000;
	}
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button) {

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button) {
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button) {

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y) {
}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h) {

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg) {

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo) {

}

