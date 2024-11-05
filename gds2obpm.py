import gdspy
import numpy as np

# GDSIIファイルの読み込み
def read_gdsii(file_path):
    gdsii = gdspy.GdsLibrary()
    gdsii.read_gds(file_path)
    return gdsii

# レイヤーごとのポリゴンを取得
def get_polygons_by_layer(gdsii):
    polygons_by_layer = {}
    for cell in gdsii.top_level():
        # get_polygons(by_spec=True)は辞書を返す
        polygons_dict = cell.get_polygons(by_spec=True)
        for (layer, datatype), polygons in polygons_dict.items():
            if layer not in polygons_by_layer:
                polygons_by_layer[layer] = []
            polygons_by_layer[layer].extend(polygons)
    return polygons_by_layer

# ポリゴンをボクセル化
def voxelize_polygons(polygons, voxel_size):
    voxel_set = set()
    for coords in polygons:
        coords = np.array(coords)
        if len(coords) == 0:
            continue
        min_x = np.min(coords[:, 0])
        max_x = np.max(coords[:, 0])
        min_y = np.min(coords[:, 1])
        max_y = np.max(coords[:, 1])
        x_range = np.arange(min_x, max_x + voxel_size, voxel_size)
        y_range = np.arange(min_y, max_y + voxel_size, voxel_size)
        for x in x_range:
            for y in y_range:
                if point_in_polygon((x, y), coords):
                    voxel_set.add((x, y))
    return voxel_set

# ポイントがポリゴン内にあるか判定（射影法）
def point_in_polygon(point, polygon):
    x, y = point
    num = len(polygon)
    j = num - 1
    odd = False
    for i in range(num):
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / ((yj - yi) + 1e-9) + xi):
            odd = not odd
        j = i
    return odd

# xmesh、ymesh、zmeshの計算
def calculate_mesh(all_voxel_data, voxel_size, z_layers):
    x_coords = [x for (x, y) in all_voxel_data]
    y_coords = [y for (x, y) in all_voxel_data]

    min_x = np.floor(min(x_coords) / voxel_size) * voxel_size
    max_x = np.ceil(max(x_coords) / voxel_size) * voxel_size
    min_y = np.floor(min(y_coords) / voxel_size) * voxel_size
    max_y = np.ceil(max(y_coords) / voxel_size) * voxel_size

    x_range = np.arange(min_x, max_x + voxel_size, voxel_size)
    y_range = np.arange(min_y, max_y + voxel_size, voxel_size)

    xmesh = f"xmesh = {min_x:.6f} {len(x_range)-1} {max_x:.6f}"
    ymesh = f"ymesh = {min_y:.6f} {len(y_range)-1} {max_y:.6f}"

    # zmeshの設定
    zmesh_values = []
    for i in range(len(z_layers)-1):
        z_start = z_layers[i]
        z_end = z_layers[i+1]
        zmesh_values.append(f"{z_start:.6f} 1 {z_end:.6f}")
    zmesh = "zmesh = " + " ".join(zmesh_values)

    return xmesh, ymesh, zmesh

# ファイル出力（指定された書式に従う）
def output_file(geometry_definitions, material_definitions, xmesh, ymesh, zmesh, output_path):
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("OpenBPM 4 2\n")
        f.write("title = GDSIIボクセル化結果\n")
        # メッシュ定義
        f.write(f"{xmesh}\n")
        f.write(f"{ymesh}\n")
        f.write(f"{zmesh}\n")
        # マテリアル定義
        for mat_id, mat_props in material_definitions.items():
            f.write(f"material = {mat_id} {mat_props['epsilon_r']} {mat_props['sigma']} {mat_props['mu_r']} {mat_props['sigma_m']}\n")
        # ジオメトリ定義
        for geom_id, geom in geometry_definitions.items():
            f.write(f"geometry = {geom['material_id']} {geom['shape']} {geom['x1']:.6f} {geom['x2']:.6f} {geom['y1']:.6f} {geom['y2']:.6f} {geom['z1']:.6f} {geom['z2']:.6f}\n")
        f.write("end\n")

# メイン処理
def main():
    gdsii_file = "input.gds"  # GDSIIファイルのパスを指定
    output_file_path = "output.txt"
    voxel_size = 0.001  # ボクセルサイズを適切に設定

    # GDSIIファイルの読み込み
    gdsii = read_gdsii(gdsii_file)

    # レイヤーごとのポリゴン取得
    polygons_by_layer = get_polygons_by_layer(gdsii)

    # マテリアル定義（例としてダミーの値を使用）
    material_definitions = {
        2: {'epsilon_r': 2.0, 'sigma': 0.0, 'mu_r': 1.0, 'sigma_m': 0.0}
    }

    # ジオメトリ定義の初期化
    geometry_definitions = {}
    geom_id = 1
    all_voxel_data = set()
    z_layers = []

    # レイヤー番号のリストを取得し、ソート
    layers = sorted(polygons_by_layer.keys())

    # レイヤーの厚さを定義（例として0.1）
    layer_thickness = 0.1

    # GDSIIレイヤー番号からZ座標へのマッピング
    layer_z_positions = {}
    current_z = 0.0
    for layer in layers:
        z_start = current_z
        z_end = current_z + layer_thickness
        layer_z_positions[layer] = (z_start, z_end)
        z_layers.append(z_start)
        z_layers.append(z_end)
        current_z = z_end

    # 各レイヤーについて処理
    for layer in layers:
        polygons = polygons_by_layer[layer]
        # ボクセル化
        voxel_data = voxelize_polygons(polygons, voxel_size)

        # 全体のボクセルデータを統合
        all_voxel_data.update(voxel_data)

        # ジオメトリ定義の作成
        if not voxel_data:
            continue  # データがない場合はスキップ

        x_coords = [x for (x, y) in voxel_data]
        y_coords = [y for (x, y) in voxel_data]

        # レイヤーのZ座標を取得
        z_start, z_end = layer_z_positions[layer]

        geometry_definitions[geom_id] = {
            'material_id': 2,  # マテリアルIDを適切に設定
            'shape': 1,  # 直方体
            'x1': min(x_coords),
            'x2': max(x_coords),
            'y1': min(y_coords),
            'y2': max(y_coords),
            'z1': z_start,
            'z2': z_end
        }
        geom_id += 1

    # Z座標の範囲を統合し、重複を削除してソート
    z_layers = sorted(set(z_layers))

    # メッシュの計算
    xmesh, ymesh, zmesh = calculate_mesh(all_voxel_data, voxel_size, z_layers)

    # ファイル出力
    output_file(geometry_definitions, material_definitions, xmesh, ymesh, zmesh, output_file_path)

if __name__ == "__main__":
    main()

