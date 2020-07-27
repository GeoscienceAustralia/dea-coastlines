<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis simplifyLocal="1" labelsEnabled="1" readOnly="0" hasScaleBasedVisibilityFlag="1" simplifyDrawingHints="0" simplifyDrawingTol="1" version="3.4.5-Madeira" simplifyMaxScale="1" maxScale="15000" simplifyAlgorithm="0" styleCategories="AllStyleCategories" minScale="500000">
  <flags>
    <Identifiable>1</Identifiable>
    <Removable>1</Removable>
    <Searchable>1</Searchable>
  </flags>
  <renderer-v2 forceraster="0" enableorderby="1" graduatedMethod="GraduatedColor" type="graduatedSymbol" symbollevels="0" attr="if(  (&quot;sig_time&quot; &lt;=0.01), rate_time, NULL)">
    <ranges>
      <range lower="-50.000000000000000" upper="-2.500000000000000" symbol="0" render="true" label=" -50.00 - -2.50 "/>
      <range lower="-2.500000000000000" upper="-1.000000000000000" symbol="1" render="true" label=" -2.50 - -1.00 "/>
      <range lower="-1.000000000000000" upper="-0.600000000000000" symbol="2" render="true" label=" -1.00 - -0.60 "/>
      <range lower="-0.600000000000000" upper="-0.300000000000000" symbol="3" render="true" label=" -0.60 - -0.30 "/>
      <range lower="-0.300000000000000" upper="-0.100000000000000" symbol="4" render="true" label=" -0.30 - -0.10 "/>
      <range lower="-0.100000000000000" upper="0.000000000000000" symbol="5" render="true" label=" -0.10 - 0.00 "/>
      <range lower="0.000000000000000" upper="0.100000000000000" symbol="6" render="true" label=" 0.00 - 0.10 "/>
      <range lower="0.100000000000000" upper="0.300000000000000" symbol="7" render="true" label=" 0.10 - 0.30 "/>
      <range lower="0.300000000000000" upper="0.600000000000000" symbol="8" render="true" label=" 0.30 - 0.60 "/>
      <range lower="0.600000000000000" upper="1.000000000000000" symbol="9" render="true" label=" 0.60 - 1.00 "/>
      <range lower="1.000000000000000" upper="2.500000000000000" symbol="10" render="true" label=" 1.00 - 2.50 "/>
      <range lower="2.500000000000000" upper="50.000000000000000" symbol="11" render="true" label=" 2.50 - 50.00 "/>
    </ranges>
    <symbols>
      <symbol type="marker" name="0" clip_to_extent="1" force_rhr="0" alpha="1">
        <layer enabled="1" locked="0" pass="0" class="SimpleMarker">
          <prop k="angle" v="0"/>
          <prop k="color" v="202,0,32,255"/>
          <prop k="horizontal_anchor_point" v="1"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="name" v="circle"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="offset_unit" v="MM"/>
          <prop k="outline_color" v="255,255,255,0"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0"/>
          <prop k="outline_width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="scale_method" v="diameter"/>
          <prop k="size" v="2.8"/>
          <prop k="size_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="size_unit" v="MM"/>
          <prop k="vertical_anchor_point" v="1"/>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
      <symbol type="marker" name="1" clip_to_extent="1" force_rhr="0" alpha="1">
        <layer enabled="1" locked="0" pass="0" class="SimpleMarker">
          <prop k="angle" v="0"/>
          <prop k="color" v="218,60,67,255"/>
          <prop k="horizontal_anchor_point" v="1"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="name" v="circle"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="offset_unit" v="MM"/>
          <prop k="outline_color" v="255,255,255,0"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0"/>
          <prop k="outline_width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="scale_method" v="diameter"/>
          <prop k="size" v="2.4"/>
          <prop k="size_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="size_unit" v="MM"/>
          <prop k="vertical_anchor_point" v="1"/>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
      <symbol type="marker" name="10" clip_to_extent="1" force_rhr="0" alpha="1">
        <layer enabled="1" locked="0" pass="0" class="SimpleMarker">
          <prop k="angle" v="0"/>
          <prop k="color" v="56,144,193,255"/>
          <prop k="horizontal_anchor_point" v="1"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="name" v="circle"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="offset_unit" v="MM"/>
          <prop k="outline_color" v="255,255,255,0"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0"/>
          <prop k="outline_width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="scale_method" v="diameter"/>
          <prop k="size" v="2.4"/>
          <prop k="size_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="size_unit" v="MM"/>
          <prop k="vertical_anchor_point" v="1"/>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
      <symbol type="marker" name="11" clip_to_extent="1" force_rhr="0" alpha="1">
        <layer enabled="1" locked="0" pass="0" class="SimpleMarker">
          <prop k="angle" v="0"/>
          <prop k="color" v="5,113,176,255"/>
          <prop k="horizontal_anchor_point" v="1"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="name" v="circle"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="offset_unit" v="MM"/>
          <prop k="outline_color" v="255,255,255,0"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0"/>
          <prop k="outline_width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="scale_method" v="diameter"/>
          <prop k="size" v="2.8"/>
          <prop k="size_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="size_unit" v="MM"/>
          <prop k="vertical_anchor_point" v="1"/>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
      <symbol type="marker" name="2" clip_to_extent="1" force_rhr="0" alpha="1">
        <layer enabled="1" locked="0" pass="0" class="SimpleMarker">
          <prop k="angle" v="0"/>
          <prop k="color" v="233,120,103,255"/>
          <prop k="horizontal_anchor_point" v="1"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="name" v="circle"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="offset_unit" v="MM"/>
          <prop k="outline_color" v="255,255,255,0"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0"/>
          <prop k="outline_width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="scale_method" v="diameter"/>
          <prop k="size" v="2"/>
          <prop k="size_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="size_unit" v="MM"/>
          <prop k="vertical_anchor_point" v="1"/>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
      <symbol type="marker" name="3" clip_to_extent="1" force_rhr="0" alpha="1">
        <layer enabled="1" locked="0" pass="0" class="SimpleMarker">
          <prop k="angle" v="0"/>
          <prop k="color" v="245,173,141,255"/>
          <prop k="horizontal_anchor_point" v="1"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="name" v="circle"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="offset_unit" v="MM"/>
          <prop k="outline_color" v="255,255,255,0"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0"/>
          <prop k="outline_width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="scale_method" v="diameter"/>
          <prop k="size" v="1.6"/>
          <prop k="size_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="size_unit" v="MM"/>
          <prop k="vertical_anchor_point" v="1"/>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
      <symbol type="marker" name="4" clip_to_extent="1" force_rhr="0" alpha="1">
        <layer enabled="1" locked="0" pass="0" class="SimpleMarker">
          <prop k="angle" v="0"/>
          <prop k="color" v="246,203,183,255"/>
          <prop k="horizontal_anchor_point" v="1"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="name" v="circle"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="offset_unit" v="MM"/>
          <prop k="outline_color" v="255,255,255,0"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0"/>
          <prop k="outline_width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="scale_method" v="diameter"/>
          <prop k="size" v="1.3"/>
          <prop k="size_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="size_unit" v="MM"/>
          <prop k="vertical_anchor_point" v="1"/>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
      <symbol type="marker" name="5" clip_to_extent="1" force_rhr="0" alpha="1">
        <layer enabled="1" locked="0" pass="0" class="SimpleMarker">
          <prop k="angle" v="0"/>
          <prop k="color" v="247,232,226,255"/>
          <prop k="horizontal_anchor_point" v="1"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="name" v="circle"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="offset_unit" v="MM"/>
          <prop k="outline_color" v="255,255,255,0"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0"/>
          <prop k="outline_width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="scale_method" v="diameter"/>
          <prop k="size" v="1"/>
          <prop k="size_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="size_unit" v="MM"/>
          <prop k="vertical_anchor_point" v="1"/>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
      <symbol type="marker" name="6" clip_to_extent="1" force_rhr="0" alpha="1">
        <layer enabled="1" locked="0" pass="0" class="SimpleMarker">
          <prop k="angle" v="0"/>
          <prop k="color" v="229,238,243,255"/>
          <prop k="horizontal_anchor_point" v="1"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="name" v="circle"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="offset_unit" v="MM"/>
          <prop k="outline_color" v="255,255,255,0"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0"/>
          <prop k="outline_width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="scale_method" v="diameter"/>
          <prop k="size" v="1"/>
          <prop k="size_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="size_unit" v="MM"/>
          <prop k="vertical_anchor_point" v="1"/>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
      <symbol type="marker" name="7" clip_to_extent="1" force_rhr="0" alpha="1">
        <layer enabled="1" locked="0" pass="0" class="SimpleMarker">
          <prop k="angle" v="0"/>
          <prop k="color" v="192,220,234,255"/>
          <prop k="horizontal_anchor_point" v="1"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="name" v="circle"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="offset_unit" v="MM"/>
          <prop k="outline_color" v="255,255,255,0"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0"/>
          <prop k="outline_width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="scale_method" v="diameter"/>
          <prop k="size" v="1.3"/>
          <prop k="size_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="size_unit" v="MM"/>
          <prop k="vertical_anchor_point" v="1"/>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
      <symbol type="marker" name="8" clip_to_extent="1" force_rhr="0" alpha="1">
        <layer enabled="1" locked="0" pass="0" class="SimpleMarker">
          <prop k="angle" v="0"/>
          <prop k="color" v="155,202,225,255"/>
          <prop k="horizontal_anchor_point" v="1"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="name" v="circle"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="offset_unit" v="MM"/>
          <prop k="outline_color" v="255,255,255,0"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0"/>
          <prop k="outline_width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="scale_method" v="diameter"/>
          <prop k="size" v="1.6"/>
          <prop k="size_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="size_unit" v="MM"/>
          <prop k="vertical_anchor_point" v="1"/>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
      <symbol type="marker" name="9" clip_to_extent="1" force_rhr="0" alpha="1">
        <layer enabled="1" locked="0" pass="0" class="SimpleMarker">
          <prop k="angle" v="0"/>
          <prop k="color" v="107,174,210,255"/>
          <prop k="horizontal_anchor_point" v="1"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="name" v="circle"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="offset_unit" v="MM"/>
          <prop k="outline_color" v="255,255,255,0"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0"/>
          <prop k="outline_width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="scale_method" v="diameter"/>
          <prop k="size" v="2"/>
          <prop k="size_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="size_unit" v="MM"/>
          <prop k="vertical_anchor_point" v="1"/>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
    </symbols>
    <source-symbol>
      <symbol type="marker" name="0" clip_to_extent="1" force_rhr="0" alpha="1">
        <layer enabled="1" locked="0" pass="0" class="SimpleMarker">
          <prop k="angle" v="0"/>
          <prop k="color" v="255,158,23,255"/>
          <prop k="horizontal_anchor_point" v="1"/>
          <prop k="joinstyle" v="bevel"/>
          <prop k="name" v="circle"/>
          <prop k="offset" v="0,0"/>
          <prop k="offset_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="offset_unit" v="MM"/>
          <prop k="outline_color" v="255,255,255,88"/>
          <prop k="outline_style" v="solid"/>
          <prop k="outline_width" v="0"/>
          <prop k="outline_width_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="outline_width_unit" v="MM"/>
          <prop k="scale_method" v="diameter"/>
          <prop k="size" v="2"/>
          <prop k="size_map_unit_scale" v="3x:0,0,0,0,0,0"/>
          <prop k="size_unit" v="MM"/>
          <prop k="vertical_anchor_point" v="1"/>
          <data_defined_properties>
            <Option type="Map">
              <Option value="" type="QString" name="name"/>
              <Option name="properties"/>
              <Option value="collection" type="QString" name="type"/>
            </Option>
          </data_defined_properties>
        </layer>
      </symbol>
    </source-symbol>
    <colorramp type="gradient" name="[source]">
      <prop k="color1" v="202,0,32,255"/>
      <prop k="color2" v="5,113,176,255"/>
      <prop k="discrete" v="0"/>
      <prop k="rampType" v="gradient"/>
      <prop k="stops" v="0.25;244,165,130,255:0.5;247,247,247,255:0.75;146,197,222,255"/>
    </colorramp>
    <mode name="equal"/>
    <symmetricMode enabled="false" symmetryPoint="0" astride="false"/>
    <rotation/>
    <sizescale/>
    <labelformat trimtrailingzeroes="false" format=" %1 - %2 " decimalplaces="2"/>
    <orderby>
      <orderByClause nullsFirst="0" asc="1"> abs("rate_time" )</orderByClause>
    </orderby>
  </renderer-v2>
  <labeling type="simple">
    <settings>
      <text-style fontLetterSpacing="0" fontFamily="Sans Serif" textColor="255,255,255,255" isExpression="1" fontCapitals="0" useSubstitutions="0" previewBkgrdColor="#ffffff" fontUnderline="0" fontSize="7.5" fontItalic="0" fontSizeMapUnitScale="3x:0,0,0,0,0,0" namedStyle="Normal" fieldName="if(abs(rate_time) > 0.1, concat(format_number(rate_time, 1), ' m (±', format_number(se_time * 1.96, 1), ')'), NULL)" fontStrikeout="0" multilineHeight="5" textOpacity="1" blendMode="0" fontSizeUnit="Point" fontWordSpacing="0" fontWeight="50">
        <text-buffer bufferDraw="1" bufferNoFill="1" bufferBlendMode="0" bufferSize="0.7778" bufferOpacity="0.5" bufferSizeUnits="MM" bufferSizeMapUnitScale="3x:0,0,0,0,0,0" bufferColor="0,0,0,255" bufferJoinStyle="128"/>
        <background shapeDraw="0" shapeBorderColor="128,128,128,0" shapeOffsetX="0" shapeRadiiMapUnitScale="3x:0,0,0,0,0,0" shapeRadiiY="5" shapeJoinStyle="64" shapeSizeType="0" shapeSizeY="5" shapeOffsetY="0" shapeOffsetUnit="MM" shapeRadiiX="5" shapeFillColor="255,255,255,0" shapeRotation="0" shapeRotationType="0" shapeSizeMapUnitScale="3x:0,0,0,0,0,0" shapeType="0" shapeBorderWidthUnit="MM" shapeOpacity="1" shapeOffsetMapUnitScale="3x:0,0,0,0,0,0" shapeBorderWidth="0" shapeBlendMode="0" shapeSizeX="5" shapeBorderWidthMapUnitScale="3x:0,0,0,0,0,0" shapeSizeUnit="MM" shapeRadiiUnit="MM" shapeSVGFile=""/>
        <shadow shadowOffsetUnit="MM" shadowOpacity="1" shadowRadiusUnit="MM" shadowUnder="0" shadowOffsetMapUnitScale="3x:0,0,0,0,0,0" shadowColor="255,255,255,255" shadowBlendMode="8" shadowDraw="1" shadowRadius="1.5" shadowScale="100" shadowOffsetDist="1" shadowRadiusAlphaOnly="0" shadowOffsetGlobal="1" shadowRadiusMapUnitScale="3x:0,0,0,0,0,0" shadowOffsetAngle="135"/>
        <substitutions/>
      </text-style>
      <text-format wrapChar="" reverseDirectionSymbol="0" placeDirectionSymbol="0" leftDirectionSymbol="&lt;" autoWrapLength="0" multilineAlign="3" formatNumbers="0" decimals="1" useMaxLineLengthForAutoWrap="1" plussign="0" addDirectionSymbol="0" rightDirectionSymbol=">"/>
      <placement repeatDistanceMapUnitScale="3x:0,0,0,0,0,0" offsetUnits="MM" preserveRotation="1" maxCurvedCharAngleOut="-25" repeatDistance="0" offsetType="0" predefinedPositionOrder="TR,TL,BR,BL,R,L,TSR,BSR" centroidInside="0" labelOffsetMapUnitScale="3x:0,0,0,0,0,0" fitInPolygonOnly="0" maxCurvedCharAngleIn="25" distMapUnitScale="3x:0,0,0,0,0,0" quadOffset="3" xOffset="-5" placementFlags="10" distUnits="MM" placement="1" repeatDistanceUnits="MM" dist="0" priority="1" centroidWhole="0" yOffset="0" rotationAngle="0"/>
      <rendering displayAll="0" limitNumLabels="0" obstacle="1" drawLabels="1" scaleMin="0" zIndex="0" maxNumLabels="50" scaleMax="1000000" scaleVisibility="1" minFeatureSize="0" labelPerPart="0" obstacleFactor="1" upsidedownLabels="0" fontMaxPixelSize="10000" mergeLines="0" fontLimitPixelSize="0" fontMinPixelSize="3" obstacleType="0"/>
      <dd_properties>
        <Option type="Map">
          <Option value="" type="QString" name="name"/>
          <Option type="Map" name="properties">
            <Option type="Map" name="Priority">
              <Option value="true" type="bool" name="active"/>
              <Option value="abs(&quot;rate_time&quot;)" type="QString" name="expression"/>
              <Option value="3" type="int" name="type"/>
            </Option>
          </Option>
          <Option value="collection" type="QString" name="type"/>
        </Option>
      </dd_properties>
    </settings>
  </labeling>
  <customproperties>
    <property value="1987-01-01" key="dualview/previewExpressions"/>
    <property value="0" key="embeddedWidgets/count"/>
    <property key="variableNames"/>
    <property key="variableValues"/>
  </customproperties>
  <blendMode>0</blendMode>
  <featureBlendMode>0</featureBlendMode>
  <layerOpacity>1</layerOpacity>
  <SingleCategoryDiagramRenderer diagramType="Histogram" attributeLegend="1">
    <DiagramCategory height="15" backgroundAlpha="255" backgroundColor="#ffffff" minScaleDenominator="15000" scaleBasedVisibility="0" sizeScale="3x:0,0,0,0,0,0" minimumSize="0" lineSizeScale="3x:0,0,0,0,0,0" scaleDependency="Area" enabled="0" maxScaleDenominator="1e+8" penAlpha="255" labelPlacementMethod="XHeight" penWidth="0" penColor="#000000" width="15" diagramOrientation="Up" barWidth="5" sizeType="MM" rotationOffset="270" opacity="1" lineSizeType="MM">
      <fontProperties description="Sans Serif,9,-1,5,50,0,0,0,0,0" style=""/>
      <attribute field="" color="#000000" label=""/>
    </DiagramCategory>
  </SingleCategoryDiagramRenderer>
  <DiagramLayerSettings linePlacementFlags="18" zIndex="0" showAll="1" placement="0" priority="0" obstacle="0" dist="0">
    <properties>
      <Option type="Map">
        <Option value="" type="QString" name="name"/>
        <Option name="properties"/>
        <Option value="collection" type="QString" name="type"/>
      </Option>
    </properties>
  </DiagramLayerSettings>
  <geometryOptions geometryPrecision="0" removeDuplicateNodes="0">
    <activeChecks/>
    <checkConfiguration/>
  </geometryOptions>
  <fieldConfiguration>
    <field name="rate_time">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="sig_time">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="se_time">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="outl_time">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_1988">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_1989">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_1990">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_1991">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_1992">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_1993">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_1994">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_1995">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_1996">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_1997">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_1998">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_1999">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2000">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2001">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2002">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2003">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2004">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2005">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2006">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2007">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2008">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2009">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2010">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2011">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2012">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2013">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2014">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2015">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2016">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2017">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2018">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="dist_2019">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="retreat">
      <editWidget type="Range">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="growth">
      <editWidget type="Range">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="sce">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="nsm">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="max_year">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
    <field name="min_year">
      <editWidget type="TextEdit">
        <config>
          <Option/>
        </config>
      </editWidget>
    </field>
  </fieldConfiguration>
  <aliases>
    <alias field="rate_time" name="" index="0"/>
    <alias field="sig_time" name="" index="1"/>
    <alias field="se_time" name="" index="2"/>
    <alias field="outl_time" name="" index="3"/>
    <alias field="dist_1988" name="" index="4"/>
    <alias field="dist_1989" name="" index="5"/>
    <alias field="dist_1990" name="" index="6"/>
    <alias field="dist_1991" name="" index="7"/>
    <alias field="dist_1992" name="" index="8"/>
    <alias field="dist_1993" name="" index="9"/>
    <alias field="dist_1994" name="" index="10"/>
    <alias field="dist_1995" name="" index="11"/>
    <alias field="dist_1996" name="" index="12"/>
    <alias field="dist_1997" name="" index="13"/>
    <alias field="dist_1998" name="" index="14"/>
    <alias field="dist_1999" name="" index="15"/>
    <alias field="dist_2000" name="" index="16"/>
    <alias field="dist_2001" name="" index="17"/>
    <alias field="dist_2002" name="" index="18"/>
    <alias field="dist_2003" name="" index="19"/>
    <alias field="dist_2004" name="" index="20"/>
    <alias field="dist_2005" name="" index="21"/>
    <alias field="dist_2006" name="" index="22"/>
    <alias field="dist_2007" name="" index="23"/>
    <alias field="dist_2008" name="" index="24"/>
    <alias field="dist_2009" name="" index="25"/>
    <alias field="dist_2010" name="" index="26"/>
    <alias field="dist_2011" name="" index="27"/>
    <alias field="dist_2012" name="" index="28"/>
    <alias field="dist_2013" name="" index="29"/>
    <alias field="dist_2014" name="" index="30"/>
    <alias field="dist_2015" name="" index="31"/>
    <alias field="dist_2016" name="" index="32"/>
    <alias field="dist_2017" name="" index="33"/>
    <alias field="dist_2018" name="" index="34"/>
    <alias field="dist_2019" name="" index="35"/>
    <alias field="retreat" name="" index="36"/>
    <alias field="growth" name="" index="37"/>
    <alias field="sce" name="" index="38"/>
    <alias field="nsm" name="" index="39"/>
    <alias field="max_year" name="" index="40"/>
    <alias field="min_year" name="" index="41"/>
  </aliases>
  <excludeAttributesWMS/>
  <excludeAttributesWFS/>
  <defaults>
    <default field="rate_time" applyOnUpdate="0" expression=""/>
    <default field="sig_time" applyOnUpdate="0" expression=""/>
    <default field="se_time" applyOnUpdate="0" expression=""/>
    <default field="outl_time" applyOnUpdate="0" expression=""/>
    <default field="dist_1988" applyOnUpdate="0" expression=""/>
    <default field="dist_1989" applyOnUpdate="0" expression=""/>
    <default field="dist_1990" applyOnUpdate="0" expression=""/>
    <default field="dist_1991" applyOnUpdate="0" expression=""/>
    <default field="dist_1992" applyOnUpdate="0" expression=""/>
    <default field="dist_1993" applyOnUpdate="0" expression=""/>
    <default field="dist_1994" applyOnUpdate="0" expression=""/>
    <default field="dist_1995" applyOnUpdate="0" expression=""/>
    <default field="dist_1996" applyOnUpdate="0" expression=""/>
    <default field="dist_1997" applyOnUpdate="0" expression=""/>
    <default field="dist_1998" applyOnUpdate="0" expression=""/>
    <default field="dist_1999" applyOnUpdate="0" expression=""/>
    <default field="dist_2000" applyOnUpdate="0" expression=""/>
    <default field="dist_2001" applyOnUpdate="0" expression=""/>
    <default field="dist_2002" applyOnUpdate="0" expression=""/>
    <default field="dist_2003" applyOnUpdate="0" expression=""/>
    <default field="dist_2004" applyOnUpdate="0" expression=""/>
    <default field="dist_2005" applyOnUpdate="0" expression=""/>
    <default field="dist_2006" applyOnUpdate="0" expression=""/>
    <default field="dist_2007" applyOnUpdate="0" expression=""/>
    <default field="dist_2008" applyOnUpdate="0" expression=""/>
    <default field="dist_2009" applyOnUpdate="0" expression=""/>
    <default field="dist_2010" applyOnUpdate="0" expression=""/>
    <default field="dist_2011" applyOnUpdate="0" expression=""/>
    <default field="dist_2012" applyOnUpdate="0" expression=""/>
    <default field="dist_2013" applyOnUpdate="0" expression=""/>
    <default field="dist_2014" applyOnUpdate="0" expression=""/>
    <default field="dist_2015" applyOnUpdate="0" expression=""/>
    <default field="dist_2016" applyOnUpdate="0" expression=""/>
    <default field="dist_2017" applyOnUpdate="0" expression=""/>
    <default field="dist_2018" applyOnUpdate="0" expression=""/>
    <default field="dist_2019" applyOnUpdate="0" expression=""/>
    <default field="retreat" applyOnUpdate="0" expression=""/>
    <default field="growth" applyOnUpdate="0" expression=""/>
    <default field="sce" applyOnUpdate="0" expression=""/>
    <default field="nsm" applyOnUpdate="0" expression=""/>
    <default field="max_year" applyOnUpdate="0" expression=""/>
    <default field="min_year" applyOnUpdate="0" expression=""/>
  </defaults>
  <constraints>
    <constraint exp_strength="0" notnull_strength="0" field="rate_time" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="sig_time" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="se_time" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="outl_time" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_1988" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_1989" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_1990" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_1991" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_1992" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_1993" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_1994" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_1995" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_1996" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_1997" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_1998" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_1999" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2000" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2001" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2002" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2003" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2004" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2005" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2006" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2007" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2008" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2009" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2010" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2011" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2012" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2013" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2014" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2015" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2016" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2017" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2018" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="dist_2019" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="retreat" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="growth" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="sce" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="nsm" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="max_year" constraints="0" unique_strength="0"/>
    <constraint exp_strength="0" notnull_strength="0" field="min_year" constraints="0" unique_strength="0"/>
  </constraints>
  <constraintExpressions>
    <constraint field="rate_time" desc="" exp=""/>
    <constraint field="sig_time" desc="" exp=""/>
    <constraint field="se_time" desc="" exp=""/>
    <constraint field="outl_time" desc="" exp=""/>
    <constraint field="dist_1988" desc="" exp=""/>
    <constraint field="dist_1989" desc="" exp=""/>
    <constraint field="dist_1990" desc="" exp=""/>
    <constraint field="dist_1991" desc="" exp=""/>
    <constraint field="dist_1992" desc="" exp=""/>
    <constraint field="dist_1993" desc="" exp=""/>
    <constraint field="dist_1994" desc="" exp=""/>
    <constraint field="dist_1995" desc="" exp=""/>
    <constraint field="dist_1996" desc="" exp=""/>
    <constraint field="dist_1997" desc="" exp=""/>
    <constraint field="dist_1998" desc="" exp=""/>
    <constraint field="dist_1999" desc="" exp=""/>
    <constraint field="dist_2000" desc="" exp=""/>
    <constraint field="dist_2001" desc="" exp=""/>
    <constraint field="dist_2002" desc="" exp=""/>
    <constraint field="dist_2003" desc="" exp=""/>
    <constraint field="dist_2004" desc="" exp=""/>
    <constraint field="dist_2005" desc="" exp=""/>
    <constraint field="dist_2006" desc="" exp=""/>
    <constraint field="dist_2007" desc="" exp=""/>
    <constraint field="dist_2008" desc="" exp=""/>
    <constraint field="dist_2009" desc="" exp=""/>
    <constraint field="dist_2010" desc="" exp=""/>
    <constraint field="dist_2011" desc="" exp=""/>
    <constraint field="dist_2012" desc="" exp=""/>
    <constraint field="dist_2013" desc="" exp=""/>
    <constraint field="dist_2014" desc="" exp=""/>
    <constraint field="dist_2015" desc="" exp=""/>
    <constraint field="dist_2016" desc="" exp=""/>
    <constraint field="dist_2017" desc="" exp=""/>
    <constraint field="dist_2018" desc="" exp=""/>
    <constraint field="dist_2019" desc="" exp=""/>
    <constraint field="retreat" desc="" exp=""/>
    <constraint field="growth" desc="" exp=""/>
    <constraint field="sce" desc="" exp=""/>
    <constraint field="nsm" desc="" exp=""/>
    <constraint field="max_year" desc="" exp=""/>
    <constraint field="min_year" desc="" exp=""/>
  </constraintExpressions>
  <expressionfields/>
  <attributeactions>
    <defaultAction value="{00000000-0000-0000-0000-000000000000}" key="Canvas"/>
  </attributeactions>
  <attributetableconfig sortExpression="&quot;move_abs&quot;" sortOrder="1" actionWidgetStyle="dropDown">
    <columns>
      <column hidden="1" type="actions" width="-1"/>
      <column hidden="0" type="field" name="rate_time" width="-1"/>
      <column hidden="0" type="field" name="sig_time" width="-1"/>
      <column hidden="0" type="field" name="outl_time" width="-1"/>
      <column hidden="0" type="field" name="se_time" width="-1"/>
      <column hidden="0" type="field" name="dist_1988" width="-1"/>
      <column hidden="0" type="field" name="dist_1989" width="-1"/>
      <column hidden="0" type="field" name="dist_1990" width="-1"/>
      <column hidden="0" type="field" name="dist_1991" width="-1"/>
      <column hidden="0" type="field" name="dist_1992" width="-1"/>
      <column hidden="0" type="field" name="dist_1993" width="-1"/>
      <column hidden="0" type="field" name="dist_1994" width="-1"/>
      <column hidden="0" type="field" name="dist_1995" width="-1"/>
      <column hidden="0" type="field" name="dist_1996" width="-1"/>
      <column hidden="0" type="field" name="dist_1997" width="-1"/>
      <column hidden="0" type="field" name="dist_1998" width="-1"/>
      <column hidden="0" type="field" name="dist_1999" width="-1"/>
      <column hidden="0" type="field" name="dist_2000" width="-1"/>
      <column hidden="0" type="field" name="dist_2001" width="-1"/>
      <column hidden="0" type="field" name="dist_2002" width="-1"/>
      <column hidden="0" type="field" name="dist_2003" width="-1"/>
      <column hidden="0" type="field" name="dist_2004" width="-1"/>
      <column hidden="0" type="field" name="dist_2005" width="-1"/>
      <column hidden="0" type="field" name="dist_2006" width="-1"/>
      <column hidden="0" type="field" name="dist_2007" width="-1"/>
      <column hidden="0" type="field" name="dist_2008" width="-1"/>
      <column hidden="0" type="field" name="dist_2009" width="-1"/>
      <column hidden="0" type="field" name="dist_2010" width="-1"/>
      <column hidden="0" type="field" name="dist_2011" width="-1"/>
      <column hidden="0" type="field" name="dist_2012" width="-1"/>
      <column hidden="0" type="field" name="dist_2013" width="-1"/>
      <column hidden="0" type="field" name="dist_2014" width="-1"/>
      <column hidden="0" type="field" name="dist_2015" width="-1"/>
      <column hidden="0" type="field" name="dist_2016" width="-1"/>
      <column hidden="0" type="field" name="dist_2017" width="-1"/>
      <column hidden="0" type="field" name="dist_2018" width="-1"/>
      <column hidden="0" type="field" name="dist_2019" width="-1"/>
      <column hidden="0" type="field" name="retreat" width="-1"/>
      <column hidden="0" type="field" name="growth" width="-1"/>
      <column hidden="0" type="field" name="sce" width="-1"/>
      <column hidden="0" type="field" name="nsm" width="-1"/>
      <column hidden="0" type="field" name="max_year" width="-1"/>
      <column hidden="0" type="field" name="min_year" width="-1"/>
    </columns>
  </attributetableconfig>
  <conditionalstyles>
    <rowstyles/>
    <fieldstyles/>
  </conditionalstyles>
  <editform tolerant="1"></editform>
  <editforminit/>
  <editforminitcodesource>0</editforminitcodesource>
  <editforminitfilepath></editforminitfilepath>
  <editforminitcode><![CDATA[# -*- coding: utf-8 -*-
"""
QGIS forms can have a Python function that is called when the form is
opened.

Use this function to add extra logic to your forms.

Enter the name of the function in the "Python Init function"
field.
An example follows:
"""
from qgis.PyQt.QtWidgets import QWidget

def my_form_open(dialog, layer, feature):
	geom = feature.geometry()
	control = dialog.findChild(QWidget, "MyLineEdit")
]]></editforminitcode>
  <featformsuppress>0</featformsuppress>
  <editorlayout>generatedlayout</editorlayout>
  <editable>
    <field name="1987-01-01" editable="1"/>
    <field name="1988" editable="1"/>
    <field name="1988-01-01" editable="1"/>
    <field name="1989" editable="1"/>
    <field name="1989-01-01" editable="1"/>
    <field name="1990" editable="1"/>
    <field name="1990-01-01" editable="1"/>
    <field name="1991" editable="1"/>
    <field name="1991-01-01" editable="1"/>
    <field name="1992" editable="1"/>
    <field name="1992-01-01" editable="1"/>
    <field name="1993" editable="1"/>
    <field name="1993-01-01" editable="1"/>
    <field name="1994" editable="1"/>
    <field name="1994-01-01" editable="1"/>
    <field name="1995" editable="1"/>
    <field name="1995-01-01" editable="1"/>
    <field name="1996" editable="1"/>
    <field name="1996-01-01" editable="1"/>
    <field name="1997" editable="1"/>
    <field name="1997-01-01" editable="1"/>
    <field name="1998" editable="1"/>
    <field name="1998-01-01" editable="1"/>
    <field name="1999" editable="1"/>
    <field name="1999-01-01" editable="1"/>
    <field name="2000" editable="1"/>
    <field name="2000-01-01" editable="1"/>
    <field name="2001" editable="1"/>
    <field name="2001-01-01" editable="1"/>
    <field name="2002" editable="1"/>
    <field name="2002-01-01" editable="1"/>
    <field name="2003" editable="1"/>
    <field name="2003-01-01" editable="1"/>
    <field name="2004" editable="1"/>
    <field name="2004-01-01" editable="1"/>
    <field name="2005" editable="1"/>
    <field name="2005-01-01" editable="1"/>
    <field name="2006" editable="1"/>
    <field name="2006-01-01" editable="1"/>
    <field name="2007" editable="1"/>
    <field name="2007-01-01" editable="1"/>
    <field name="2008" editable="1"/>
    <field name="2008-01-01" editable="1"/>
    <field name="2009" editable="1"/>
    <field name="2009-01-01" editable="1"/>
    <field name="2010" editable="1"/>
    <field name="2010-01-01" editable="1"/>
    <field name="2011" editable="1"/>
    <field name="2011-01-01" editable="1"/>
    <field name="2012" editable="1"/>
    <field name="2012-01-01" editable="1"/>
    <field name="2013" editable="1"/>
    <field name="2013-01-01" editable="1"/>
    <field name="2014" editable="1"/>
    <field name="2014-01-01" editable="1"/>
    <field name="2015" editable="1"/>
    <field name="2015-01-01" editable="1"/>
    <field name="2016" editable="1"/>
    <field name="2016-01-01" editable="1"/>
    <field name="2017" editable="1"/>
    <field name="2017-01-01" editable="1"/>
    <field name="2018" editable="1"/>
    <field name="2018-01-01" editable="1"/>
    <field name="breakpoints" editable="1"/>
    <field name="dist_1988" editable="1"/>
    <field name="dist_1989" editable="1"/>
    <field name="dist_1990" editable="1"/>
    <field name="dist_1991" editable="1"/>
    <field name="dist_1992" editable="1"/>
    <field name="dist_1993" editable="1"/>
    <field name="dist_1994" editable="1"/>
    <field name="dist_1995" editable="1"/>
    <field name="dist_1996" editable="1"/>
    <field name="dist_1997" editable="1"/>
    <field name="dist_1998" editable="1"/>
    <field name="dist_1999" editable="1"/>
    <field name="dist_2000" editable="1"/>
    <field name="dist_2001" editable="1"/>
    <field name="dist_2002" editable="1"/>
    <field name="dist_2003" editable="1"/>
    <field name="dist_2004" editable="1"/>
    <field name="dist_2005" editable="1"/>
    <field name="dist_2006" editable="1"/>
    <field name="dist_2007" editable="1"/>
    <field name="dist_2008" editable="1"/>
    <field name="dist_2009" editable="1"/>
    <field name="dist_2010" editable="1"/>
    <field name="dist_2011" editable="1"/>
    <field name="dist_2012" editable="1"/>
    <field name="dist_2013" editable="1"/>
    <field name="dist_2014" editable="1"/>
    <field name="dist_2015" editable="1"/>
    <field name="dist_2016" editable="1"/>
    <field name="dist_2017" editable="1"/>
    <field name="dist_2018" editable="1"/>
    <field name="dist_2019" editable="1"/>
    <field name="eln_outl" editable="1"/>
    <field name="eln_rate" editable="1"/>
    <field name="eln_sig" editable="1"/>
    <field name="growth" editable="1"/>
    <field name="incpt_IOD" editable="1"/>
    <field name="incpt_IPO" editable="1"/>
    <field name="incpt_PDO" editable="1"/>
    <field name="incpt_SAM" editable="1"/>
    <field name="incpt_SOI" editable="1"/>
    <field name="incpt_tide" editable="1"/>
    <field name="incpt_time" editable="1"/>
    <field name="lan_outl" editable="1"/>
    <field name="lan_rate" editable="1"/>
    <field name="lan_sig" editable="1"/>
    <field name="layer" editable="1"/>
    <field name="max_year" editable="1"/>
    <field name="min_year" editable="1"/>
    <field name="mov_outl" editable="1"/>
    <field name="mov_rate" editable="1"/>
    <field name="mov_sig" editable="1"/>
    <field name="move_abs" editable="1"/>
    <field name="neg_outl" editable="1"/>
    <field name="neg_rate" editable="1"/>
    <field name="neg_sig" editable="1"/>
    <field name="nsm" editable="1"/>
    <field name="orig_ogc_fid" editable="1"/>
    <field name="outl_IOD" editable="1"/>
    <field name="outl_IPO" editable="1"/>
    <field name="outl_PDO" editable="1"/>
    <field name="outl_SAM" editable="1"/>
    <field name="outl_SOI" editable="1"/>
    <field name="outl_tide" editable="1"/>
    <field name="outl_time" editable="1"/>
    <field name="path" editable="1"/>
    <field name="pos_outl" editable="1"/>
    <field name="pos_rate" editable="1"/>
    <field name="pos_sig" editable="1"/>
    <field name="rate_IOD" editable="1"/>
    <field name="rate_IPO" editable="1"/>
    <field name="rate_PDO" editable="1"/>
    <field name="rate_SAM" editable="1"/>
    <field name="rate_SOI" editable="1"/>
    <field name="rate_tide" editable="1"/>
    <field name="rate_time" editable="1"/>
    <field name="retreat" editable="1"/>
    <field name="sce" editable="1"/>
    <field name="se_time" editable="1"/>
    <field name="sig_IOD" editable="1"/>
    <field name="sig_IPO" editable="1"/>
    <field name="sig_PDO" editable="1"/>
    <field name="sig_SAM" editable="1"/>
    <field name="sig_SOI" editable="1"/>
    <field name="sig_tide" editable="1"/>
    <field name="sig_time" editable="1"/>
    <field name="soi_outl" editable="1"/>
    <field name="soi_rate" editable="1"/>
    <field name="soi_sig" editable="1"/>
  </editable>
  <labelOnTop>
    <field name="1987-01-01" labelOnTop="0"/>
    <field name="1988" labelOnTop="0"/>
    <field name="1988-01-01" labelOnTop="0"/>
    <field name="1989" labelOnTop="0"/>
    <field name="1989-01-01" labelOnTop="0"/>
    <field name="1990" labelOnTop="0"/>
    <field name="1990-01-01" labelOnTop="0"/>
    <field name="1991" labelOnTop="0"/>
    <field name="1991-01-01" labelOnTop="0"/>
    <field name="1992" labelOnTop="0"/>
    <field name="1992-01-01" labelOnTop="0"/>
    <field name="1993" labelOnTop="0"/>
    <field name="1993-01-01" labelOnTop="0"/>
    <field name="1994" labelOnTop="0"/>
    <field name="1994-01-01" labelOnTop="0"/>
    <field name="1995" labelOnTop="0"/>
    <field name="1995-01-01" labelOnTop="0"/>
    <field name="1996" labelOnTop="0"/>
    <field name="1996-01-01" labelOnTop="0"/>
    <field name="1997" labelOnTop="0"/>
    <field name="1997-01-01" labelOnTop="0"/>
    <field name="1998" labelOnTop="0"/>
    <field name="1998-01-01" labelOnTop="0"/>
    <field name="1999" labelOnTop="0"/>
    <field name="1999-01-01" labelOnTop="0"/>
    <field name="2000" labelOnTop="0"/>
    <field name="2000-01-01" labelOnTop="0"/>
    <field name="2001" labelOnTop="0"/>
    <field name="2001-01-01" labelOnTop="0"/>
    <field name="2002" labelOnTop="0"/>
    <field name="2002-01-01" labelOnTop="0"/>
    <field name="2003" labelOnTop="0"/>
    <field name="2003-01-01" labelOnTop="0"/>
    <field name="2004" labelOnTop="0"/>
    <field name="2004-01-01" labelOnTop="0"/>
    <field name="2005" labelOnTop="0"/>
    <field name="2005-01-01" labelOnTop="0"/>
    <field name="2006" labelOnTop="0"/>
    <field name="2006-01-01" labelOnTop="0"/>
    <field name="2007" labelOnTop="0"/>
    <field name="2007-01-01" labelOnTop="0"/>
    <field name="2008" labelOnTop="0"/>
    <field name="2008-01-01" labelOnTop="0"/>
    <field name="2009" labelOnTop="0"/>
    <field name="2009-01-01" labelOnTop="0"/>
    <field name="2010" labelOnTop="0"/>
    <field name="2010-01-01" labelOnTop="0"/>
    <field name="2011" labelOnTop="0"/>
    <field name="2011-01-01" labelOnTop="0"/>
    <field name="2012" labelOnTop="0"/>
    <field name="2012-01-01" labelOnTop="0"/>
    <field name="2013" labelOnTop="0"/>
    <field name="2013-01-01" labelOnTop="0"/>
    <field name="2014" labelOnTop="0"/>
    <field name="2014-01-01" labelOnTop="0"/>
    <field name="2015" labelOnTop="0"/>
    <field name="2015-01-01" labelOnTop="0"/>
    <field name="2016" labelOnTop="0"/>
    <field name="2016-01-01" labelOnTop="0"/>
    <field name="2017" labelOnTop="0"/>
    <field name="2017-01-01" labelOnTop="0"/>
    <field name="2018" labelOnTop="0"/>
    <field name="2018-01-01" labelOnTop="0"/>
    <field name="breakpoints" labelOnTop="0"/>
    <field name="dist_1988" labelOnTop="0"/>
    <field name="dist_1989" labelOnTop="0"/>
    <field name="dist_1990" labelOnTop="0"/>
    <field name="dist_1991" labelOnTop="0"/>
    <field name="dist_1992" labelOnTop="0"/>
    <field name="dist_1993" labelOnTop="0"/>
    <field name="dist_1994" labelOnTop="0"/>
    <field name="dist_1995" labelOnTop="0"/>
    <field name="dist_1996" labelOnTop="0"/>
    <field name="dist_1997" labelOnTop="0"/>
    <field name="dist_1998" labelOnTop="0"/>
    <field name="dist_1999" labelOnTop="0"/>
    <field name="dist_2000" labelOnTop="0"/>
    <field name="dist_2001" labelOnTop="0"/>
    <field name="dist_2002" labelOnTop="0"/>
    <field name="dist_2003" labelOnTop="0"/>
    <field name="dist_2004" labelOnTop="0"/>
    <field name="dist_2005" labelOnTop="0"/>
    <field name="dist_2006" labelOnTop="0"/>
    <field name="dist_2007" labelOnTop="0"/>
    <field name="dist_2008" labelOnTop="0"/>
    <field name="dist_2009" labelOnTop="0"/>
    <field name="dist_2010" labelOnTop="0"/>
    <field name="dist_2011" labelOnTop="0"/>
    <field name="dist_2012" labelOnTop="0"/>
    <field name="dist_2013" labelOnTop="0"/>
    <field name="dist_2014" labelOnTop="0"/>
    <field name="dist_2015" labelOnTop="0"/>
    <field name="dist_2016" labelOnTop="0"/>
    <field name="dist_2017" labelOnTop="0"/>
    <field name="dist_2018" labelOnTop="0"/>
    <field name="dist_2019" labelOnTop="0"/>
    <field name="eln_outl" labelOnTop="0"/>
    <field name="eln_rate" labelOnTop="0"/>
    <field name="eln_sig" labelOnTop="0"/>
    <field name="growth" labelOnTop="0"/>
    <field name="incpt_IOD" labelOnTop="0"/>
    <field name="incpt_IPO" labelOnTop="0"/>
    <field name="incpt_PDO" labelOnTop="0"/>
    <field name="incpt_SAM" labelOnTop="0"/>
    <field name="incpt_SOI" labelOnTop="0"/>
    <field name="incpt_tide" labelOnTop="0"/>
    <field name="incpt_time" labelOnTop="0"/>
    <field name="lan_outl" labelOnTop="0"/>
    <field name="lan_rate" labelOnTop="0"/>
    <field name="lan_sig" labelOnTop="0"/>
    <field name="layer" labelOnTop="0"/>
    <field name="max_year" labelOnTop="0"/>
    <field name="min_year" labelOnTop="0"/>
    <field name="mov_outl" labelOnTop="0"/>
    <field name="mov_rate" labelOnTop="0"/>
    <field name="mov_sig" labelOnTop="0"/>
    <field name="move_abs" labelOnTop="0"/>
    <field name="neg_outl" labelOnTop="0"/>
    <field name="neg_rate" labelOnTop="0"/>
    <field name="neg_sig" labelOnTop="0"/>
    <field name="nsm" labelOnTop="0"/>
    <field name="orig_ogc_fid" labelOnTop="0"/>
    <field name="outl_IOD" labelOnTop="0"/>
    <field name="outl_IPO" labelOnTop="0"/>
    <field name="outl_PDO" labelOnTop="0"/>
    <field name="outl_SAM" labelOnTop="0"/>
    <field name="outl_SOI" labelOnTop="0"/>
    <field name="outl_tide" labelOnTop="0"/>
    <field name="outl_time" labelOnTop="0"/>
    <field name="path" labelOnTop="0"/>
    <field name="pos_outl" labelOnTop="0"/>
    <field name="pos_rate" labelOnTop="0"/>
    <field name="pos_sig" labelOnTop="0"/>
    <field name="rate_IOD" labelOnTop="0"/>
    <field name="rate_IPO" labelOnTop="0"/>
    <field name="rate_PDO" labelOnTop="0"/>
    <field name="rate_SAM" labelOnTop="0"/>
    <field name="rate_SOI" labelOnTop="0"/>
    <field name="rate_tide" labelOnTop="0"/>
    <field name="rate_time" labelOnTop="0"/>
    <field name="retreat" labelOnTop="0"/>
    <field name="sce" labelOnTop="0"/>
    <field name="se_time" labelOnTop="0"/>
    <field name="sig_IOD" labelOnTop="0"/>
    <field name="sig_IPO" labelOnTop="0"/>
    <field name="sig_PDO" labelOnTop="0"/>
    <field name="sig_SAM" labelOnTop="0"/>
    <field name="sig_SOI" labelOnTop="0"/>
    <field name="sig_tide" labelOnTop="0"/>
    <field name="sig_time" labelOnTop="0"/>
    <field name="soi_outl" labelOnTop="0"/>
    <field name="soi_rate" labelOnTop="0"/>
    <field name="soi_sig" labelOnTop="0"/>
  </labelOnTop>
  <widgets/>
  <previewExpression>1987-01-01</previewExpression>
  <mapTip></mapTip>
  <layerGeometryType>0</layerGeometryType>
</qgis>
