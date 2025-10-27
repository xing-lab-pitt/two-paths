import numpy as np
import pandas as pd
import os
import tifffile as tiff
import pandas as pd
import plotly.graph_objects as go

class Config:
    
    # main directory
    home_dir = ''

    # phate representation file path
    phate_representation_path = f'{home_dir}/output/phate_representation.csv'

    # used for training the manifold
    feature_names = [
                     'cyclina_mean_ratio', 
                     'cyclinb_mean_ratio', 
                     'cyclind_mean_ratio', 
                     'dapi_nucleus_total_intensity', 
                     'vimentin_cell_total_intensity'
                     ]  
    # feature for clustering
    clustering_feature = 'dapi_nucleus_total_intensity'

    # assign group id, only for multiple experiments
    group_id_dict = {'14days_treated': 0, '3days_treated': 1, '3days_no_treated': 2}
    
    # save plotted figures under output/figure/out_folder_name
    out_folder_name = 'self_selected'   


    # extra info record cell identity, if multiple experiments, add group
    extra_info = ['cell_id', 'tile_id', 'nucleus_id', 'bbox']  

    # replace z axis feature name
    replace_z = 'vimentin_cell_total_intensity'   # when umap components is 2
    # replace_z = None  # when umap components is 3



def plot_group(merged_df, feature,  group_id_dict, save_dir):

    inverted_group_id_dict = dict(zip(group_id_dict.values(), group_id_dict.keys()))

    fig_list = []

    color_map = ['viridis', 'Cividis', 'Plasma']

    color_map_index = 0

    for group_value in list(group_id_dict.values()):
        
        data_df = merged_df[merged_df['group'] == group_value]
        cell_id = data_df['cell_id'].to_numpy().astype(np.uint16)
        tile_id = data_df['tile_id'].to_numpy()
        nucleus_id = data_df['nucleus_id'].to_numpy().astype(np.uint16)
        bbox = data_df['bbox'].to_numpy()

        # use 3 sigma on colorbar
        channel_values = data_df[feature].to_numpy()
        
        mean_value = np.mean(channel_values)
        std_dev = np.std(channel_values)

        # calculate range from -3σ to +3σ 
        lower_bound = mean_value - 3 * std_dev
        upper_bound = mean_value + 3 * std_dev

        # 3d scatter
        fig = go.Scatter3d(
            x=data_df['x'],
            y=data_df['y'],
            z=data_df['z'],
            mode='markers',
            name = inverted_group_id_dict[group_value],
            marker=dict(
                size=5, 
                color=channel_values,
                colorscale=color_map[color_map_index],
                colorbar=dict(
                    title=feature, 
                    x=0.8,      
                    tickmode='array',
                    tickvals=[lower_bound, mean_value, upper_bound],  
                    ticktext=[f'{np.min(channel_values):.2f}', f'{np.mean(channel_values):.2f}', f'{np.max(channel_values):.2f}']    # use sigma value to color but show the original value
                ),
                cmin=lower_bound,  # minvalue of colorbar
                cmax=upper_bound,  # maxvalue 0f colorbar
                opacity=0.8,
                showscale=True,
            ),
            customdata=[
                [f'Cell id: {cell_id[i]}', 
                f'Tile id: {tile_id[i]}', 
                f'Cell_nucleus_id: {nucleus_id[i]}', 
                f'Bbox: {bbox[i]}', 
                f'Value: {channel_values[i]}',
                ] 
                for i in range(len(cell_id))
            ],
            hovertemplate="""
            <b>Custom Info</b><br>
            X: %{x}<br>
            Y: %{y}<br>
            Z: %{z}<br>
            %{customdata[0]}<br>
            %{customdata[1]}<br>
            %{customdata[2]}<br>
            %{customdata[3]}<br>
            %{customdata[4]}<br>
            <extra></extra>
            """,
            hoverlabel=dict(
                font=dict(
                    size=10
                )
            )
        )

        fig_list.append(fig)

        color_map_index += 1
        
    fig = go.Figure(data = fig_list)

    fig.update_layout( 
        title='3D PhatePlot',
        scene=dict(
            xaxis_title='Phate1',
            yaxis_title='Phate2',
            zaxis_title='Phate3' if Config.replace_z is None else Config.replace_z
        ),
        margin=dict(l=0, r=0, b=0, t=0),
        legend=dict(x=1, y=0.9, font=dict(size=12), itemsizing='constant')  
        )


    # Embed the additional HTML, add colorbar
    custom_html = """
    <script>
    document.addEventListener('DOMContentLoaded', function(){

        let legendClicked = false;

        function updateColorbar(plot) {
            setTimeout(function(){
                let visibleCount = 0;
                let visibleIndex = -1;
                // check visible state of each trace
                for (let i = 0; i < plot.data.length; i++) {
                    if (plot.data[i].visible !== 'legendonly' && plot.data[i].visible !== false) {
                        visibleCount++;
                        visibleIndex = i;
                    }
                }
             
                let newShowscale = [false, false, false];
                if (visibleCount === 1) {
                    newShowscale[visibleIndex] = true;
                }
                Plotly.restyle(plot, {'marker.showscale': newShowscale}, [0, 1, 2]);

            }, 100); // delay
        }

        var plots = document.getElementsByClassName('js-plotly-plot');

        // Attach to hover event
        Array.from(plots).forEach(function(plot) {

            plot.on('plotly_legendclick', function() {
                legendClicked = true;
            });

            plot.on('plotly_afterplot', function() {
                if (legendClicked) {
                    updateColorbar(plot);
                    legendClicked = false;
                }
            });
        });
    });
    </script>

    """

    # Write the figure to an HTML file
    with open(save_dir, "w") as f:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))  # Export the figure
        f.write(custom_html)  # Append the custom HTML and JavaScript



def plot_clustering(merged_df, feature,  group_id_dict, save_dir):

    inverted_group_id_dict = dict(zip(group_id_dict.values(), group_id_dict.keys()))
    fig_list = []

    color_map = [[[0, '#EB929B'], [1, '#43a8a0']],   
                 [[0, '#EB929B'], [1, '#43a8a0']],
                 [[0, '#EB929B'], [1, '#43a8a0']],]


    color_map_index = 0

    for group_value in list(group_id_dict.values()):
        
        data_df = merged_df[merged_df['group'] == group_value]
        cell_id = data_df['cell_id'].to_numpy().astype(np.uint16)
        tile_id = data_df['tile_id'].to_numpy()
        nucleus_id = data_df['nucleus_id'].to_numpy().astype(np.uint16)
        bbox = data_df['bbox'].to_numpy()

        channel_values = data_df[feature].to_numpy()
        
        # clustering on channel values 2 clusters
        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters=2, random_state=0).fit(channel_values.reshape(-1, 1))
        channel_values = kmeans.labels_
        # label count check and make sure label 1 is the larger one
        if np.sum(channel_values == 1) < np.sum(channel_values == 0):
            channel_values = 1 - channel_values

        # 3d scatter
        fig = go.Scatter3d(
            x=data_df['x'],
            y=data_df['y'],
            z=data_df['z'],
            mode='markers',
            name = inverted_group_id_dict[group_value],
            marker=dict(
                size=5, 
                color=channel_values,
                colorscale=color_map[color_map_index],
                opacity=0.8,
                showscale=False,
            ),  
            customdata=[
                [f'Cell id: {cell_id[i]}', 
                f'Tile id: {tile_id[i]}', 
                f'Cell_nucleus_id: {nucleus_id[i]}', 
                f'Bbox: {bbox[i]}', 
                f'Value: {channel_values[i]}',
                ] 
                for i in range(len(cell_id))
            ],
            hovertemplate="""
            <b>Custom Info</b><br>
            X: %{x}<br>
            Y: %{y}<br>
            Z: %{z}<br>
            %{customdata[0]}<br>
            %{customdata[1]}<br>
            %{customdata[2]}<br>
            %{customdata[3]}<br>
            %{customdata[4]}<br>
            <extra></extra>
            """,
            hoverlabel=dict(
                font=dict(
                    size=10
                )
            )
        )
            

        fig_list.append(fig)

        color_map_index += 1
        
    fig = go.Figure(data = fig_list)

    fig.update_layout( 
        title='3D PhatePlot',
        scene=dict(
            xaxis_title='Phate1',
            yaxis_title='Phate2',
            zaxis_title='Phate3' if Config.replace_z is None else Config.replace_z
        ),
        margin=dict(l=0, r=0, b=0, t=0),
        legend=dict(x=1, y=0.9, font=dict(size=12), itemsizing='constant')  
        )

    fig.write_html(save_dir)



# plot all experiments
def plot_group_umap_on_selected_features(feature_names): 

    merged_df = pd.read_csv(f'/home/zhz187/Nikon/combined_three/output/phate_representation.csv')
    
    for feature in feature_names:
        if os.path.exists(f'{Config.home_dir}/output/figure/{Config.out_folder_name}/{feature}.html'):
            continue
        
        plot_group(merged_df, feature, Config.group_id_dict, save_dir = f'{Config.home_dir}/output/figure/{Config.out_folder_name}/{feature}.html')

    if not os.path.exists(f'{Config.home_dir}/output/figure/{Config.out_folder_name}/{Config.clustering_feature}_clustering.html'):
        plot_clustering(merged_df, Config.clustering_feature, Config.group_id_dict, save_dir = f'{Config.home_dir}/output/figure/{Config.out_folder_name}/{Config.clustering_feature}_clustering.html')


if __name__ == '__main__':
    
    plot_group_umap_on_selected_features(Config.feature_names)
