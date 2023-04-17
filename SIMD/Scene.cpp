#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    Vector3f dir = { 0.0,0.0,0.0 };
    Vector3f indir = { 0.0,0.0,0.0 };
    //1.�ж��Ƿ��н��㣺�����볡���������ཻ
    Intersection inter = Scene::intersect(ray);
    //���û����
    if (!inter.happened) {
        return dir;//return 0,0,0
    }
    //2.ray�򵽹�Դ�ˣ�˵����Ⱦ����ֻ����ǰ����Է�������ֱ�ӷ��ز��ʵ��Է�����
    if (inter.m->hasEmission()) {
        if (depth == 0) {//��һ�δ򵽹�
            return inter.m->getEmission();
        }
        else return dir;//����򵽹⣬ֱ�ӷ���0��0.0

    }
    //3.ray�����壺���ʱ��ſ�ʼ����α�������Ĳ���

    //�Գ����еĹ�Դ���в������õ�������light_pos��pdf_light
    Intersection light_pos;
    float pdf_light = 0.0f;
    sampleLight(light_pos, pdf_light);

    //3.1����ֱ�ӹ���

    //�����һЩ����
    Vector3f p = inter.coords;
    Vector3f N = inter.normal.normalized();
    Vector3f wo = ray.direction;//����ָ�򳡾�
    //��Դ��һЩ����
    Vector3f xx = light_pos.coords;
    Vector3f NN = light_pos.normal.normalized();
    Vector3f ws = (p - xx).normalized();//��Դָ������
    float dis = (p - xx).norm();//���߾���
    float dis2 = dotProduct((p - xx), (p - xx));

    //�жϹ�Դ��������Ƿ����ڵ���
    //����һ�����ߣ�����Ϊws ��Դxx -> ����p
    Ray light_to_obj(xx, ws);//Ray(orig,dir)
    Intersection light_to_scene = Scene::intersect(light_to_obj);
    //����dis>light_to_scene.distance��˵�����ڵ�����ô���Ÿ��������ɣ�
    if (light_to_scene.happened && (light_to_scene.distance - dis > -EPSILON)) {//û���ڵ�
        //Ϊ�˸�����α���룬���趨һЩ����
        Vector3f L_i = light_pos.emit;//��ǿ
        Vector3f f_r = inter.m->eval(wo, -ws, N);//���ʣ�����˵�ˣ�BRDF==���ʣ�ws���������
        float cos_theta = dotProduct(-ws, N);//����н�
        float cos_theta_l = dotProduct(ws, NN);//��Դ�н�
        dir = L_i * f_r * cos_theta * cos_theta_l / dis2 / pdf_light;
    }

    //3.2��ӹ���

    //����˹���̶�
    //Scene.hpp���Ѿ�������P_RR:RussianRoulette=0.8
    float ksi = get_random_float();//���ȡ[0,1]
    if (ksi < RussianRoulette) {
        //�����ӹ���

        //�������һ��wi����
        Vector3f wi = inter.m->sample(wo, N).normalized();//�����wi��ʵû������㣬���ص���һ������ķ���
        Ray r(p, wi);
        Intersection obj_to_scene = Scene::intersect(r);
        //����������&&���岻�ǹ�Դ
        if (obj_to_scene.happened && !obj_to_scene.m->hasEmission()) {
            Vector3f f_r = inter.m->eval(wo, wi, N);//wo���������
            float cos_theta = dotProduct(wi, N);
            float pdf_hemi = inter.m->pdf(wo, wi, N);
            indir = castRay(r, depth + 1) * f_r * cos_theta / pdf_hemi / RussianRoulette;
        }
    }
    return dir + indir;
}